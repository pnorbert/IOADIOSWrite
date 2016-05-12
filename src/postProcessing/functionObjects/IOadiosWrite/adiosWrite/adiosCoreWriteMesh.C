/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
     \\/     M anipulation  |
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

#include "adiosCoreWrite.H"
#include "polyMesh.H"


// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //

void Foam::adiosCoreWrite::definePatchAttributes(const polyMesh& mesh)
{
    const polyPatchList& patches = mesh.boundaryMesh();
    stringList pNames(patches.size());
    stringList pTypes(patches.size());
    forAll(patches, patchI)
    {
        const polyPatch& p = patches[patchI];

        pNames[patchI] = p.name();
        pTypes[patchI] = p.type();
    }

    const fileName varPath = adiosCore::regionPath(mesh.name());

    defineIntAttribute("nPatches",     varPath, patches.size());
    defineListAttribute("patch-names", varPath, pNames);
    defineListAttribute("patch-types", varPath, pTypes);
}


void Foam::adiosCoreWrite::defineMeshPoints(const polyMesh& mesh)
{
    const fileName varPath = adiosCore::meshPath(mesh.name());

    // polyMesh/nPoints (summary)
    // polyMesh/points: 2D array (N points x 3 coordinates)
    defineIntVariable(varPath/"nPoints");
    defineVectorVariable(varPath/"points",  mesh.nPoints());
}


void Foam::adiosCoreWrite::defineMeshFaces(const polyMesh& mesh)
{
    const fileName varPath = adiosCore::meshPath(mesh.name());

    defineIntVariable(varPath/"nCells");        // polyMesh/nCells
    defineIntVariable(varPath/"nFaces");        // polyMesh/nFaces
    defineIntVariable(varPath/"nInternalFaces"); // polyMesh/nInternalFaces

    // polyMesh/faces - save in compact form
    // need transcription for output
    {
        const faceList& faces = mesh.faces();

        // indices = nFaces+1
        label count = faces.size()+1;

        defineIntVariable(varPath/"faces"/"indices", count);

        // count size for compact format
        count = 0;
        forAll(faces, faceI)
        {
            count += faces[faceI].size();
        }

        defineIntVariable(varPath/"faces"/"content", count);
    }

    // polyMesh/owner - direct write
    defineIntVariable(varPath/"owner", mesh.faceOwner().size());

    // polyMesh/neighbour - direct write
    defineIntVariable(varPath/"neighbour", mesh.faceNeighbour().size());

    // polyMesh/boundary - byte-stream
    {
        OutputCounter os(adiosCore::strFormat);
        os << mesh.boundaryMesh();

        defineStreamVariable(varPath/"boundary", os.size());
    }
}


void Foam::adiosCoreWrite::defineMesh(const polyMesh& mesh)
{
    defineMeshPoints(mesh);
    defineMeshFaces(mesh);
}


void Foam::adiosCoreWrite::writeMeshPoints(const polyMesh& mesh)
{
    const fileName varPath = adiosCore::meshPath(mesh.name());

    writeIntVariable(varPath/"nPoints", mesh.nPoints());
    writeVariable(varPath/"points",     mesh.points());
}


void Foam::adiosCoreWrite::writeMeshFaces(const polyMesh& mesh)
{
    const fileName varPath = adiosCore::meshPath(mesh.name());

    writeIntVariable(varPath/"nCells",  mesh.nCells());
    writeIntVariable(varPath/"nFaces",  mesh.nFaces());
    writeIntVariable(varPath/"nInternalFaces", mesh.nInternalFaces());

    // use iobuffer_ to avoid reallocations
    // Needs rework?
    iobuffer_.reserve(adiosCoreWrite::maxSize());

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


void Foam::adiosCoreWrite::writeMesh(const polyMesh& mesh)
{
    writeMeshPoints(mesh);
    writeMeshFaces(mesh);
}


// ************************************************************************* //
