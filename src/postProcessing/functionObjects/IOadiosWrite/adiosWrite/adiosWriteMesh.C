/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |               2015 Norbert Podhorszki
                            |               2016 OpenCFD Ltd
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
#include "cellModeller.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::cellModel* Foam::adiosWrite::unknownModel =
    Foam::cellModeller::lookup("unknown");

const Foam::cellModel* Foam::adiosWrite::tetModel =
    Foam::cellModeller::lookup("tet");

const Foam::cellModel* Foam::adiosWrite::pyrModel =
    Foam::cellModeller::lookup("pyr");

const Foam::cellModel* Foam::adiosWrite::prismModel =
    Foam::cellModeller::lookup("prism");

const Foam::cellModel* Foam::adiosWrite::hexModel =
    Foam::cellModeller::lookup("hex");

const Foam::cellModel* Foam::adiosWrite::wedgeModel =
    Foam::cellModeller::lookup("wedge");


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

size_t Foam::adiosWrite::meshDefine(regionInfo& r)
{
    OCompactCountStream os(adiosCore::strFormat);
    size_t maxLen = 0;

    Info<< "adiosWrite::meshDefine: region"
        << r.index_ << "=" << r.name_ << " at time " << time_.timeName() << endl;

    const fvMesh& mesh = time_.lookupObject<fvMesh>(r.name_);


    // Find over all (global) number of cells per process
    r.nCells_[Pstream::myProcNo()] = mesh.cells().size();
    Pstream::gatherList(r.nCells_);
    Pstream::scatterList(r.nCells_);

    fileName varPath = "mesh" + Foam::name(r.index_);
    Info<< "varpath = " << varPath << endl;

    defineVariable(varPath/"time",      adios_integer);
    defineVariable(varPath/"timeidx",   adios_integer);

    varPath = "region" + Foam::name(r.index_) / "polyMesh";
    // Define mesh

    // polyMesh/points
    meshDefinePoints(mesh, r);

    // polyMesh/faces - save in compact form
    {
        const faceList& faces = mesh.faces();
        size_t bufLen = 0;

        // indices = nFaces+1
        label count = faces.size()+1;

        bufLen = defineVariable
        (
            varPath/"faces"/"indices",
            adios_integer,
            count
        );

        // will need transcription for output
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
            adios_integer,
            count
        );
        maxLen = Foam::max(maxLen, bufLen);
    }

    // polyMesh/owner
    defineVariable
    (
        varPath/"owner",
        adios_integer,
        mesh.faceOwner().size()
    );

    // polyMesh/neighbour
    defineVariable
    (
        varPath/"neighbour",
        adios_integer,
        mesh.faceNeighbour().size()
    );

    // polyMesh/boundary - byte-stream
    {
        os.rewind();
        os << mesh.boundaryMesh();

        size_t bufLen = os.size();
        maxLen = Foam::max(maxLen, bufLen);

        defineVariable(varPath/"boundary", adios_unsigned_byte, bufLen);
    }

    meshDefineCells(mesh, r);

    //meshDefineBoundaries();

//
// constant/polyMesh/pointProcAddressing
// constant/polyMesh/boundaryProcAddressing
// constant/polyMesh/cellProcAddressing
// constant/polyMesh/faceProcAddressing

#if 1
    Info<< "Writing polyMesh" << endl;
    const_cast<fvMesh&>(mesh).setInstance(mesh.time().timeName());
    mesh.polyMesh::write();
#endif
    // cells -> compactIOList

// faceCompactIOList faces_;
// typedef CompactIOList<face, label> faceCompactIOList;
//
// list of labels (starting points)
// list of elements of the basetype (eg, labels)

    Info<< endl;

    return maxLen;
}


void Foam::adiosWrite::meshWrite(const regionInfo& r)
{
    Info<< "adiosWrite::meshWrite:" << endl;

    const fvMesh& mesh = time_.lookupObject<fvMesh>(r.name_);

    fileName varPath = "mesh" + Foam::name(r.index_);
    Info<< "varpath = " << varPath << endl;

    int t = mesh.time().timeOutputValue();
    writeVariable(varPath/"time", &t);

    t = mesh.time().timeIndex();
    writeVariable(varPath/"timeidx", &t);

    // Write mesh
    meshWritePoints(mesh, r);
    meshWriteCells(mesh, r);

    varPath = "region" + Foam::name(r.index_) / "polyMesh";

    // polyMesh/faces - save in compact form
#if 1
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
        Pout<< "indices " << start.size() << endl;

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

        writeVariable(varPath/"faces"/"indices", start.cdata());

        // this is not particularly elegant, but should be stable
        List<label> elems(start[start.size()-1]);

        // use iobuffer_ to avoid reallocations
        // UList<label> elems
        // (
        //     reinterpret_cast<label*>(iobuffer_.data()),
        //     iobuffer_.capacity() / sizeof(label)
        // );

        Pout<< "size: " << elems.size() << endl;

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

        writeVariable(varPath/"faces"/"content", elems.cdata());
    }
#endif

    // polyMesh/owner
    writeVariable(varPath/"owner", mesh.faceOwner().cdata());

    // polyMesh/neighbour
    writeVariable(varPath/"neighbour", mesh.faceNeighbour().cdata());

    // polyMesh/boundary
    {
        OCompactBufStream os(iobuffer_, adiosCore::strFormat);
        os << mesh.boundaryMesh();

        writeVariable(varPath/"boundary", iobuffer_.cdata());
    }

    Info<< endl;
}

// face-io (compact form)
//
// sizes:
//    start-list: nFace+1
//    input-elem: ???
//
//         // Convert to compact format
//         labelList start(L.size()+1);
//
//         start[0] = 0;
//         for (label i = 1; i < start.size(); i++)
//         {
//             label prev = start[i-1];
//             start[i] = prev+L[i-1].size();
//
//             if (start[i] < prev)
//             {
//                 FatalIOErrorInFunction(os)
//                     << "Overall number of elements " << start[i]
//                     << " of CompactIOList of size "
//                     << L.size() << " overflows the representation of a label"
//                     << endl << "Please recompile with a larger representation"
//                     << " for label" << exit(FatalIOError);
//             }
//         }
//
//         List<BaseType> elems(start[start.size()-1]);
//
//         label elemI = 0;
//         forAll(L, i)
//         {
//             const T& subList = L[i];
//
//             forAll(subList, j)
//             {
//                 elems[elemI++] = subList[j];
//             }
//         }
//         os << start << elems;
//     }



size_t Foam::adiosWrite::meshDefinePoints(const fvMesh& m, regionInfo& r)
{
    Info<< "  meshDefinePoints" << endl;

    const pointField& points = m.points();

    // Find out how many points each process has
///     List<label> nPoints(Pstream::nProcs());
///     nPoints[Pstream::myProcNo()] = points.size();
///     Pstream::gatherList(nPoints);
///     Pstream::scatterList(nPoints);

    fileName varPath("mesh" + Foam::name(r.index_));
    Info<< "varpath = " << varPath << endl;

    // Define a 1D array to store number of points on each processor
    const string global = Foam::name(Pstream::nProcs());      // form a global 1D array of this info
    const string offset = Foam::name(Pstream::myProcNo());    // offsets of this process in the 1D array
    adios_define_var
    (
        groupID_,
        (varPath/"npoints").c_str(),
        "",
        adios_integer,
        "1",
        global.c_str(),
        offset.c_str()
    );
    outputSize_ += 4; // size of adios_integer

    Info<< "  " << varPath/"npoints" << "global: " << global << "offset:" << offset << endl;

    // local dimension of the 2D sub-array (n points x 3 coordinates)
    defineVectorVariable(varPath/"points", ADIOS_SCALAR, points.size());

    varPath = "region" + Foam::name(r.index_) / "polyMesh";

    // local dimension of the 2D sub-array (N points x 3 coordinates)
    defineVectorVariable(varPath/"points", ADIOS_SCALAR, points.size());

    return 0;
}


size_t Foam::adiosWrite::meshDefineCells(const fvMesh& m, regionInfo& r)
{
    Info<< "  meshDefineCells" << endl;

    // Map shapes OpenFOAM->XDMF
    Map<label> shapeLookupIndex;
    shapeLookupIndex.insert(hexModel->index(), 9);
    shapeLookupIndex.insert(prismModel->index(), 8);
    shapeLookupIndex.insert(pyrModel->index(), 7);
    shapeLookupIndex.insert(tetModel->index(), 6);
    shapeLookupIndex.insert(wedgeModel->index(), 5);
    shapeLookupIndex.insert(unknownModel->index(), 0);

    const cellList& cells  = m.cells();
    const cellShapeList& shapes = m.cellShapes();


    // Find dataset length for this process
    // this will possible give a little overhead w.r.t. storage, but on a
    // hex-dominated mesh, this is OK.
    int j = 0;
    forAll(cells, cellId)
    {
        const cellShape& shape = shapes[cellId];
        label mapIndex = shape.model().index();

        Info<< "  cell ID = " << cellId << " index = " << mapIndex
            << ". Shape name = " << shape.model().name()
            << ". Shape nPoints = " << shape.model().nPoints()
            << ". Shape nEdges = " << shape.model().nEdges()
            << ". Shape nFaces = " << shape.model().nFaces()
            << endl;

        // A registered primitive type
        if (shapeLookupIndex.found(mapIndex))
        {
            const labelList& vrtList = shapes[cellId];
            j++; // one element for shapeID;
            forAll(vrtList, i)
            {
                j++; // one element for each vertex of the cell (vrtList[i])
            }
        }

        // If the cell is empty, leave it alone
        // FIXME: I have no idea if it's okay to ignore it
        //else if  (shape.model().nPoints()==0) {
        //    Info<< "  empty cell ignored, id=" << cellId << endl;
        //}
        // If the cell is not a basic type, exit with an error
        else
        {
            //j++; // we will add a 0 type and no vertices to continue
            //continue;
            FatalErrorIn
            (
                "adiosWrite::meshWriteCells()"
            )   << "Unsupported or unknown cell type for cell number "
                << cellId
                << endl
                << exit(FatalError);
        }
    }


    // We don't need this: Find out how many points each process has
    r.cellDataSizes_[Pstream::myProcNo()] = j;
    //Pstream::gatherList(r.cellDataSizes_);
    //Pstream::scatterList(r.cellDataSizes_);


    // Create the different datasets (needs to be done collectively)
    char datasetName[80];
    char ldimstr[16];
    char gdimstr[16];
    char offstr[16];

    // Define a 1D array to store number of cells of each processor
    //sprintf (datasetName, "MESH/%s/nCells", m.time().timeName().c_str());
    sprintf (datasetName, "mesh%d/ncells", r.index_);
    sprintf (gdimstr, "%d", Pstream::nProcs());      // form a global 1D array of this info
    sprintf (offstr,  "%d", Pstream::myProcNo());      // offsets of this process in the 1D array
    adios_define_var(groupID_, datasetName, "", adios_integer, "1", gdimstr, offstr);
    outputSize_ += 4; // size of adios_integer

    // Define a local 1D array for the Cell dataset for this process
    /*sprintf
        (
            datasetName,
            "MESH/%s/processor%i/CELLS",
            m.time().timeName().c_str(),
            Pstream::myProcNo()
        );
    */
    sprintf(datasetName, "mesh%d/cells", r.index_);
    sprintf(ldimstr, "%d", r.cellDataSizes_[Pstream::myProcNo()]);
    adios_define_var(groupID_, datasetName, "", adios_integer, ldimstr, ldimstr, "0");
    outputSize_ += r.cellDataSizes_[Pstream::myProcNo()] * 4; // size of adios_integer

    return 0;
}


size_t Foam::adiosWrite::meshDefineBoundaries(const fvMesh& m, regionInfo& r)
{
    /*
    Info<< "  meshDefineBoundaries" << endl;
    const polyPatchList& patches = m.boundaryMesh();

    Info<< "----------------" << nl
        << "Patches" << nl
        << "----------------" << nl;

    forAll(patches, patchI)
    {
        const polyPatch& p = patches[patchI];
        const polyBoundaryMesh& b = p.boundaryMesh();
        const polyMesh& m = b.mesh();
        const cellList& cells  = m.cells();
        const cellShapeList& shapes = m.cellShapes();

        Info<< "  " << "patch " << patchI
            << " (start: " << p.start()
            << " size: " << p.size()
            << ") name: " << p.name()
            << nl;
    }
    */

    return 0;
}


void Foam::adiosWrite::meshWritePoints(const fvMesh& m, const regionInfo& r)
{
    Info<< "  meshWritePoints" << endl;

    fileName varPath = "mesh" + Foam::name(r.index_);

    const pointField& points = m.points();

    // Create a simple array of points (to pass on to adios_write)
    ioScalar pointList[points.size()][3];
    forAll(points, ptI)
    {
        pointList[ptI][0] = points[ptI].x();
        pointList[ptI][1] = points[ptI].y();
        pointList[ptI][2] = points[ptI].z();
    }

    // Find out how many points each process has
    List<label> nPoints(Pstream::nProcs());
    nPoints[Pstream::myProcNo()] = points.size();
    Pstream::gatherList(nPoints);
    Pstream::scatterList(nPoints);

    int n = points.size();
    writeVariable(varPath/"npoints", &n);
    Pout<<"npoints: " << n << endl;

    // Select correct dataset for this process
    /*sprintf
        (
            datasetName,
            "MESH/%s/processor%i/POINTS",
            m.time().timeName().c_str(),
            Pstream::myProcNo()
        );*/
    writeVariable(varPath/"points", pointList);

    varPath = "region" + Foam::name(r.index_) / "polyMesh";

    Info<<"writing with " << adios_type_to_string(ADIOS_IOSCALAR_TYPE) << endl;

#ifdef ADIOS_USE_SINGLE
    Info<<"using transcribe for points" << endl;
    writeVariable(varPath/"points", pointList);
#else
    writeVariable(varPath/"points", points.cdata());
#endif
}


void Foam::adiosWrite::meshWriteCells(const fvMesh& m, const regionInfo& r)
{
    Info<< "  meshWriteCells" << endl;

    // Map shapes OpenFOAM->XDMF
    Map<label> shapeLookupIndex;
    shapeLookupIndex.insert(hexModel->index(), 9);
    shapeLookupIndex.insert(prismModel->index(), 8);
    shapeLookupIndex.insert(pyrModel->index(), 7);
    shapeLookupIndex.insert(tetModel->index(), 6);
    shapeLookupIndex.insert(wedgeModel->index(), 5);
    shapeLookupIndex.insert(unknownModel->index(), 0);

    const cellList& cells  = m.cells();
    const cellShapeList& shapes = m.cellShapes();

    // Find dataset length for this process and fill dataset in one operation
    // this will possible give a little overhead w.r.t. storage, but on a
    // hex-dominated mesh, this is OK.
    int j = 0;
    int myDataset[9*cells.size()];
    forAll(cells, cellId)
    {
        const cellShape& shape = shapes[cellId];
        label mapIndex = shape.model().index();

        // A registered primitive type
        if (shapeLookupIndex.found(mapIndex))
        {
            label shapeId = shapeLookupIndex[mapIndex];
            const labelList& vrtList = shapes[cellId];

            myDataset[j] = shapeId; j++;
            forAll(vrtList, i)
            {
                myDataset[j] = vrtList[i]; j++;
            }
        }

        // If the cell is not a basic type, exit with an error
        else
        {
            //myDataset[j] = 0; j++;
            //continue;

            FatalErrorInFunction
                << "Unsupported or unknown cell type for cell number "
                << cellId << endl
                << exit(FatalError);

        }
    }


    fileName varPath("mesh" + Foam::name(r.index_));

    fileName datasetName("mesh" + Foam::name(r.index_)/"ncells");
    int s = r.cellDataSizes_[Pstream::myProcNo()];
    writeVariable(varPath/"ncells", &s);

    // Write a separate 1D array for the Cell dataset for this process
    /* sprintf(
            datasetName,
            "MESH/%s/processor%i/CELLS",
            m.time().timeName().c_str(),
            Pstream::myProcNo()
        );*/
    writeVariable(varPath/"cells", myDataset);
}


// ************************************************************************* //
