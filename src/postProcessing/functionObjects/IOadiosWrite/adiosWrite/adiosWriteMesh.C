 /*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |               2015 Norbert Podhorszki
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::adiosWrite::meshDefine(label regionID)
{
    Info<< "adiosWrite::meshDefine: region " << regions_[regionID].name_ << endl;
    
    regionInfo& r = regions_[regionID];
    const fvMesh& m = time_.lookupObject<fvMesh>(r.name_);

    // Find over all (global) number of cells per process
    r.nCells_[Pstream::myProcNo()] = m.cells().size();
    Pstream::gatherList(r.nCells_);
    Pstream::scatterList(r.nCells_);

    // Define time and time index as single scalar (written from rank 0)
    char datasetName[80];
    sprintf (datasetName, "mesh%d/time", regionID);
    adios_define_var (groupID_, datasetName,    "", adios_integer, NULL, NULL, NULL);
    sprintf (datasetName, "mesh%d/timeidx", regionID);
    adios_define_var (groupID_, datasetName, "", adios_integer, NULL, NULL, NULL);
    outputSize_ += 2*4; // size of adios_integer

    // Define mesh
    meshDefinePoints(m, regionID);
    meshDefineCells(m, regionID);
    //meshDefineBoundaries();
    
    Info<< endl;
}

void Foam::adiosWrite::meshWrite(label regionID)
{
    Info<< "adiosWrite::meshWrite:" << endl;
    
    const fvMesh& m = time_.lookupObject<fvMesh>(regions_[regionID].name_);

    char datasetName[80];
    sprintf (datasetName, "mesh%d/time", regionID);
    int t = m.time().timeOutputValue();
    adios_write(fileID_, datasetName, &t);

    sprintf (datasetName, "mesh%d/timeidx", regionID);
    t = m.time().timeIndex();
    adios_write(fileID_, datasetName, &t);

    // Write mesh
    meshWritePoints(m, regionID);
    meshWriteCells(m, regionID);
    

    Info<< endl;
}

void Foam::adiosWrite::meshDefinePoints(const fvMesh& m, label regionID)
{
    Info<< "  meshDefinePoints" << endl;

    const pointField& points = m.points();

    // Find out how many points each process has
    List<label> nPoints(Pstream::nProcs());
    nPoints[Pstream::myProcNo()] = points.size();
    Pstream::gatherList(nPoints);
    Pstream::scatterList(nPoints);

    /*
    // Sum total number of particles on all processes
    label nTotalPoints = sum(nPoints);

    int myoffset = 0;
    for (label proc=0; proc < Pstream::myProcNo(); proc++)
    {
        myoffset += nPoints[proc];
    }
    */

    // Create the dataset for points
    char datasetName[80];
    char ldimstr[16];
    char gdimstr[16];
    char offstr[16];

    // Define a 1D array to store number of points on each processor
    sprintf (datasetName, "mesh%d/npoints", regionID);
    sprintf (gdimstr, "%d", Pstream::nProcs());      // form a global 1D array of this info
    sprintf (offstr,  "%d", Pstream::myProcNo());    // offsets of this process in the 1D array
    adios_define_var (groupID_, datasetName, "", adios_integer, "1", gdimstr, offstr);
    outputSize_ += 4; // size of adios_integer

    /*
    sprintf
    (
    		datasetName,
    		"MESH/%s/processor%i/POINTS",
    		m.time().timeName().c_str(),
    		Pstream::myProcNo()
    );
    */
    sprintf (datasetName, "mesh%d/points", regionID);
    sprintf (ldimstr, "%d,3", points.size()); // local dimension of the 2D sub-array (n points x 3 coordinates)
    //sprintf (gdimstr, "%d,3", nTotalPoints);  // global dimension of the 2D array (N points x 3 coordinates)
    //sprintf (offstr,  "%d,0", myoffset);      // offsets of this process in the 2D array
    //adios_define_var (groupID_, datasetName, "", ADIOS_SCALAR, ldimstr, gdimstr, offstr);
    adios_define_var (groupID_, datasetName, "", ADIOS_SCALAR, ldimstr, "", "");

    // count the total size we are going to write from this process
    outputSize_ += points.size() * 3 * sizeof(ioScalar);
}


void Foam::adiosWrite::meshDefineCells(const fvMesh& m, label regionID)
{
    Info<< "  meshDefineCells" << endl;

    regionInfo& r = regions_[regionID];

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
    sprintf (datasetName, "mesh%d/ncells", regionID);
    sprintf (gdimstr, "%d", Pstream::nProcs());      // form a global 1D array of this info
    sprintf (offstr,  "%d", Pstream::myProcNo());      // offsets of this process in the 1D array
    adios_define_var (groupID_, datasetName, "", adios_integer, "1", gdimstr, offstr);
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
    sprintf (datasetName, "mesh%d/cells", regionID);
    sprintf (ldimstr, "%d", r.cellDataSizes_[Pstream::myProcNo()]);
    adios_define_var (groupID_, datasetName, "", adios_integer, ldimstr, ldimstr, "0");
    outputSize_ += r.cellDataSizes_[Pstream::myProcNo()] * 4; // size of adios_integer
 }

void Foam::adiosWrite::meshDefineBoundaries(const fvMesh& m, label regionID)
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
}


void Foam::adiosWrite::meshWritePoints(const fvMesh& m, label regionID)
{   
    Info<< "  meshWritePoints" << endl;
    
    const pointField& points = m.points();
    char datasetName[80];
    
    // Create a simple array of points (to pass on to adios_write)
    ioScalar pointList[points.size()][3];
    forAll(points, ptI)
    {
        pointList[ptI][0] = points[ptI].x();
        pointList[ptI][1] = points[ptI].y();
        pointList[ptI][2] = points[ptI].z();
    }
    
    sprintf (datasetName, "mesh%d/npoints", regionID);
    int n = points.size();
    adios_write(fileID_, datasetName, &n);

    // Select correct dataset for this process
    /*sprintf
        (
            datasetName,
            "MESH/%s/processor%i/POINTS",
            m.time().timeName().c_str(),
            Pstream::myProcNo()
        );*/
    sprintf (datasetName, "mesh%d/points", regionID);
    adios_write(fileID_, datasetName, pointList);
}


void Foam::adiosWrite::meshWriteCells(const fvMesh& m, label regionID)
{
    Info<< "  meshWriteCells" << endl;

    regionInfo& r = regions_[regionID];

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
            
            FatalErrorIn
            (
                "adiosWrite::meshWriteCells()"
            )   << "Unsupported or unknown cell type for cell number "
                << cellId << endl
                << exit(FatalError);
            
        }
    }
    
    
    char datasetName[80];

    sprintf (datasetName, "mesh%d/ncells", regionID);
    int s =  r.cellDataSizes_[Pstream::myProcNo()];
    adios_write (fileID_, datasetName, &s);

    // Write a separate 1D array for the Cell dataset for this process
    /*sprintf
        (
            datasetName,
            "MESH/%s/processor%i/CELLS",
            m.time().timeName().c_str(),
            Pstream::myProcNo()
        );*/
    sprintf (datasetName, "mesh%d/cells", regionID);
    adios_write (fileID_, datasetName, myDataset);
}


const Foam::cellModel* Foam::adiosWrite::unknownModel = Foam::cellModeller::
lookup
(
    "unknown"
);


const Foam::cellModel* Foam::adiosWrite::tetModel = Foam::cellModeller::
lookup
(
    "tet"
);


const Foam::cellModel* Foam::adiosWrite::pyrModel = Foam::cellModeller::
lookup
(
    "pyr"
);


const Foam::cellModel* Foam::adiosWrite::prismModel = Foam::cellModeller::
lookup
(
    "prism"
);


const Foam::cellModel* Foam::adiosWrite::hexModel = Foam::cellModeller::
lookup
(
    "hex"
);

const Foam::cellModel* Foam::adiosWrite::wedgeModel = Foam::cellModeller::
lookup
(
    "wedge"
);

// ************************************************************************* //
