/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 Norbert Podhorszki
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "adiosWrite.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

size_t Foam::adiosWrite::meshDefine(regionInfo& r)
{
    OutputCounter os(adiosCore::strFormat);
    size_t bufLen = 0;
    size_t maxLen = 0;

    Info<< "adiosWrite::meshDefine: region"
        << r.index_ << "=" << r.name_ << " at time "
        << time_.timeName() << endl;

    const fvMesh& mesh = time_.lookupObject<fvMesh>(r.name_);

#ifdef FOAM_ADIOS_CELL_SHAPES
    meshDefineCellShapes(mesh, r);
#endif

    fileName varPath = "region" + Foam::name(r.index_) / "polyMesh";

    // summary information
    defineIntVariable(varPath/"nPoints");       // polyMesh/nPoints
    defineIntVariable(varPath/"nCells");        // polyMesh/nCells
    defineIntVariable(varPath/"nFaces");        // polyMesh/nFaces
    defineIntVariable(varPath/"nInternalFaces"); // polyMesh/nInternalFaces

    // polyMesh/points: 2D array (N points x 3 coordinates)
    defineVectorVariable
    (
        varPath/"points",
        adiosTraits<scalar>::adiosType,
        mesh.nPoints()
    );

    // polyMesh/faces - save in compact form
    {
        const faceList& faces = mesh.faces();

        // indices = nFaces+1
        label count = faces.size()+1;

        bufLen = defineVariable
        (
            varPath/"faces"/"indices",
            adiosTraits<label>::adiosType,
            count
        );

        // need transcription for output
        maxLen = Foam::max(maxLen, bufLen);

        // count size for compact format
        count = 0;
        forAll(faces, faceI)
        {
            count += faces[faceI].size();
        }

        bufLen = defineVariable
        (
            varPath/"faces"/"content",
            adiosTraits<label>::adiosType,
            count
        );
        maxLen = Foam::max(maxLen, bufLen);
    }

    // polyMesh/owner - direct write
    defineVariable
    (
        varPath/"owner",
        adiosTraits<label>::adiosType,
        mesh.faceOwner().size()
    );

    // polyMesh/neighbour - direct write
    defineVariable
    (
        varPath/"neighbour",
        adiosTraits<label>::adiosType,
        mesh.faceNeighbour().size()
    );

    // polyMesh/boundary - byte-stream
    {
        os.rewind();
        os << mesh.boundaryMesh();
        bufLen = os.size();

        defineVariable(varPath/"boundary", adios_unsigned_byte, bufLen);
        maxLen = Foam::max(maxLen, bufLen);
    }

    return maxLen;
}


void Foam::adiosWrite::meshWrite(const regionInfo& r)
{
    Info<< "adiosWrite::meshWrite: region"
        << r.index_ << "=" << r.name_ << " at time "
        << time_.timeName() << endl;

    const fvMesh& mesh = time_.lookupObject<fvMesh>(r.name_);

#ifdef FOAM_ADIOS_CELL_SHAPES
    meshWriteCellShapes(mesh, r);
#endif

    fileName varPath = "region" + Foam::name(r.index_) / "polyMesh";

    writeIntVariable(varPath/"nPoints", mesh.nPoints());
    writeIntVariable(varPath/"nCells",  mesh.nCells());
    writeIntVariable(varPath/"nFaces",  mesh.nFaces());
    writeIntVariable(varPath/"nInternalFaces", mesh.nInternalFaces());

    // polyMesh/points
    writeVariable(varPath/"points", mesh.points());

    // polyMesh/faces - save in compact form
    // for this we use two separate lists
    {
        const faceList& faces = mesh.faces();
        const label nFaces = faces.size();

        // this is not particularly elegant, but should be stable
        List<label> start(nFaces+1);

        // // use iobuffer_ to avoid reallocations
        // UList<label> start
        // (
        //     reinterpret_cast<label*>(iobuffer_.data()),
        //     iobuffer_.capacity() / sizeof(label)
        // );

        // use iobuffer_ to avoid reallocations
        // Pout<< "indices " << start.size() << endl;

        // calculate start as per CompactIOList.C:
        start[0] = 0;
        forAll(faces, faceI)
        {
            label prev = start[faceI];
            start[faceI+1] = prev + faces[faceI].size();

            if (start[faceI+1] < prev)
            {
                FatalErrorInFunction
                    << "Overall number of elements " << start[faceI+1]
                    << " of CompactIOList of size "
                    << nFaces << " overflows the representation of a label"
                    << endl << "Please recompile with a larger representation"
                    << " for label" << exit(FatalIOError);
            }
        }

        writeVariable(varPath/"faces"/"indices", start);

        // this is not particularly elegant, but should be stable
        List<label> elems(start[start.size()-1]);

        // use iobuffer_ to avoid reallocations
        // UList<label> elems
        // (
        //     reinterpret_cast<label*>(iobuffer_.data()),
        //     iobuffer_.capacity() / sizeof(label)
        // );

        // Pout<< "size: " << elems.size() << endl;

        // calculate content as per CompactIOList.C:
        label elemI = 0;
        forAll(faces, faceI)
        {
            const face& f = faces[faceI];

            forAll(f, i)
            {
                elems[elemI++] = f[i];
            }
        }

        writeVariable(varPath/"faces"/"content", elems);
    }

    // polyMesh/owner
    writeVariable(varPath/"owner", mesh.faceOwner());

    // polyMesh/neighbour
    writeVariable(varPath/"neighbour", mesh.faceNeighbour());

    // polyMesh/boundary - as bytes to avoid issues with adios_string
    {
        OutputBufStreamer os(iobuffer_, adiosCore::strFormat);
        os << mesh.boundaryMesh();

        writeVariable(varPath/"boundary", iobuffer_);
    }
}


// ************************************************************************* //
