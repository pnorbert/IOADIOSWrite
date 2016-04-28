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
#include "cellModeller.H"

// * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * * //

// file-local
static inline Foam::fileName oldMeshPath(Foam::label index)
{
    return "mesh" + Foam::name(index);
}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Map<Foam::label>& Foam::adiosWrite::shapeLookupMap() const
{
    if (!shapeLookupPtr_)
    {
        shapeLookupPtr_ = new Map<label>(10);

        shapeLookupPtr_->insert
        (
            cellModeller::lookup("hex")->index(),
            9
        );
        shapeLookupPtr_->insert
        (
            cellModeller::lookup("prism")->index(),
            8
        );
        shapeLookupPtr_->insert
        (
            cellModeller::lookup("pyr")->index(),
            7
        );
        shapeLookupPtr_->insert
        (
            cellModeller::lookup("tet")->index(),
            6
        );
        shapeLookupPtr_->insert
        (
            cellModeller::lookup("wedge")->index(),
            5
        );
        shapeLookupPtr_->insert
        (
            cellModeller::lookup("unknown")->index(),
            0
        );
    }

    return *shapeLookupPtr_;
}


#ifdef FOAM_ADIOS_CELL_SHAPES

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::adiosWrite::meshDefineCellShapes(const fvMesh& mesh, regionInfo& r)
{
    Info<< "  meshDefineCellShapes" << endl;

    fileName varPath = oldMeshPath(r.index_);

    // defineIntVariable(varPath/"timeidx");
    // defineScalarVariable(varPath/"time");

    // 1D array with number of points on each processor
    defineIntVariable(varPath/"nPoints");

    // 1D array with number of cells on each processor
    defineIntVariable(varPath/"nCells");

    // polyMesh/points: 2D array (N points x 3 coordinates)
    defineVectorVariable(varPath/"points", mesh.nPoints());


    // Map shapes OpenFOAM->XDMF
    const Map<label>& shapeLookupIndex = shapeLookupMap();
    const cellShapeList& shapes = mesh.cellShapes();

    // Find dataset length for this process
    int j = 0;
    forAll(shapes, cellId)
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

            ++j;                 // shape-id;
            j += vrtList.size(); // each cell vertex
        }
        else
        {
            // If the cell is empty, leave it alone
            // FIXME: I have no idea if it's okay to ignore it
            //else if  (shape.model().nPoints()==0) {
            //    Info<< "  empty cell ignored, id=" << cellId << endl;
            //}
            // If the cell is not a basic type, exit with an error

            //j++; // we will add a 0 type and no vertices to continue
            //continue;
            FatalErrorInFunction
                << "Unsupported or unknown cell type for cell number "
                << cellId
                << endl
                << exit(FatalError);
        }
    }

    r.cellDataSizes_[Pstream::myProcNo()] = j;

    // We don't need this: Find out how many points each process has
    //Pstream::gatherList(r.cellDataSizes_);
    //Pstream::scatterList(r.cellDataSizes_);

    defineIntVariable(varPath/"cells", r.cellDataSizes_[Pstream::myProcNo()]);
}


void Foam::adiosWrite::meshWriteCellShapes(const fvMesh& mesh, const regionInfo& r)
{
    Info<< "  meshWriteCellShapes" << endl;

    fileName varPath = oldMeshPath(r.index_);

    // writeIntVariable(varPath/"timeidx", mesh.time().timeIndex());
    // writeScalarVariable(varPath/"time", mesh.time().timeOutputValue());

    // Write mesh
    writeIntVariable(varPath/"nPoints", mesh.nPoints());
    writeIntVariable(varPath/"nCells",  mesh.nCells());
    writeVariable(varPath/"points",     mesh.points());

    // Map shapes OpenFOAM->XDMF
    const Map<label>& shapeLookupIndex = shapeLookupMap();
    const cellShapeList& shapes = mesh.cellShapes();

    int j = 0;
    label labelBuffer[r.cellDataSizes_[Pstream::myProcNo()]];
    forAll(shapes, cellId)
    {
        const cellShape& shape = shapes[cellId];
        label mapIndex = shape.model().index();

        // A registered primitive type
        if (shapeLookupIndex.found(mapIndex))
        {
            labelBuffer[j++] = shapeLookupIndex[mapIndex];

            const labelList& vrtList = shapes[cellId];
            forAll(vrtList, i)
            {
                labelBuffer[j++] = vrtList[i];
            }
        }
        else
        {
            // If the cell is not a basic type, exit with an error
            // labelBuffer[j++] = 0;
            // continue;

            FatalErrorInFunction
                << "Unsupported or unknown cell type for cell number "
                << cellId << endl
                << exit(FatalError);
        }
    }

    writeVariable(varPath/"cells", labelBuffer);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif /* FOAM_ADIOS_CELL_SHAPES */


// ************************************************************************* //
