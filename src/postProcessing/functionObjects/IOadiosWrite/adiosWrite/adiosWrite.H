/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |               2015 Norbert Podhorszki
                            |               2016 OpenCFD Ltd.
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

Class
    Foam::adiosWrite

Description
    Writes selected (vector and scalar) fields and particle clouds to an
    ADIOS file.
    Based on the IOh5Write module created by
       Håkon Strandenes, hakostra@stud.ntnu.no

SourceFiles
    adiosWrite.C
    IOadiosWrite.H

Author
    Norbert Podhorszki, pnorbert@ornl.gov, 2015

\*---------------------------------------------------------------------------*/

#ifndef adiosWrite_H
#define adiosWrite_H

#include "adiosCore.H"

#include "runTimeSelectionTables.H"
#include "interpolation.H"
#include "fvMesh.H"
#include "surfaceMesh.H"
#include "fvsPatchFields.H"

#include "adios.h"
#include "adios_read.h"
#include <stdio.h> // sprintf

// If OpenFOAM is compiled with single precision, the IO will be single
// precision.
#if defined(WM_DP)
#   define ADIOS_SCALAR adios_double
    typedef double ioScalar;
#elif defined(WRITE_SP)
#   define ADIOS_SCALAR adios_float
    typedef float ioScalar;
#else
#   error "Neither WM_DP nor WM_SP is defined when compiling adiosWrite.H"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
//class fvMesh;

/*---------------------------------------------------------------------------*\
                         Class adiosWrite Declaration
\*---------------------------------------------------------------------------*/

class adiosWrite
:
    public adiosCore
{
protected:

    // Protected classes

        //- Class template used for grouping field types
        template<class Type>
        class fieldGroup
        :
            public DynamicList<word>
        {
        public:
            //- Construct null
            fieldGroup()
            :
                DynamicList<word>(0)
            {}

        };

        // Per-region info variables grouped together to be able to create a list of them
        class regionInfo
        {
        public:
            regionInfo() {} ;
            ~regionInfo() {};

            word name_;

            //- Names of volume fields per region to write
            wordReList objectNames_;

            //- Names of clouds to write
            wordReList cloudNames_;

            //- Which attributes to write for particle data
            wordReList cloudAttribs_;

            //- Categorized scalar/vector/tensor fields
            fieldGroup<scalar> scalarFields_;
            fieldGroup<vector> vectorFields_;
            fieldGroup<scalar> surfaceScalarFields_;

            //- Array containing number of cells per process (often used list)
            List<label> nCells_;

            //- Number of data points in output Cell array on each process
            //  calculated in meshDefineCells, used in meshWriteCells
            //  only the local size is used currently
            List<label> cellDataSizes_;

            //- Array containing number of particles per process (often used list)
            List<label> nParticles_;

            //- Total number of particles (sum of nParticles_ list)
            label nTotalParticles_;
        };


    // Protected data

        // Reference to the object registry
        const objectRegistry& obr_;

        //- Reference to the primary mesh database
        const fvMesh& primaryMesh_; //  primaryMesh_(refCast<const fvMesh>(obr))

        //- Reference to the time database
        const Time& time_;  // = primaryMesh_.time()

        //- List of regions holding all data needed (mesh, names)
        List<regionInfo> regions_;

        //- ADIOS I/O method to use
        word adiosMethod_;

        //- ADIOS I/O method parameters
        string methodParams_;

        //- Write interval (in number of iterations)
        label writeInterval_;

        // Restart time requested
        scalar restartTime_;

        // Used to track when to write
        label nextWrite_;
        label timeSteps_;

        // ADIOS File and other identifiers
        int64_t  fileID_;     // file handler from open()
        int64_t  groupID_;    // group of output variables
        uint64_t outputSize_; // number of bytes we write in one timestep by the process


        // MPI parallel code variables
        MPI_Comm comm_; // Pstream:: does not provide communicator
        // rank = Pstream::myProcNo()
        // number of processors = Pstream::nProcs()
        // process id: Pstream::procID (int procNo), 1 for rank 0


    // Protected Member Functions

        // General functions

            // Create and open dataset for all data
            void open();

            // Close dataset
            void close();


        // Functions for handling field data

            //- Field define wrapper
            //  Return the maximum associated stream length
            size_t fieldDefine(label regionID);

            //- Field define for a given type
            //  Return the maximum associated stream length
            template<class FieldType>
            size_t fieldDefine
            (
                const fvMesh& mesh,
                const fieldGroup<typename FieldType::value_type>& fields,
                const label regionID
            );

            //- Field write wrapper
            void fieldWrite(label regionID);

            //- Write data for a given type
            template<class FieldType>
            void fieldWrite
            (
                const fvMesh& mesh,
                const fieldGroup<typename FieldType::value_type>& fields,
                const label regionID
            );


      // Functions for handling particle cloud data

            //- Cloud define wrapper
            void cloudDefine(label regionID);

            //- Cloud write wrapper
            void cloudWrite(label regionID);

            //- Write attributes for cloud particles
            void cloudWriteAttrib
            (
                label myParticles,
                label offset,
                label nTot,
                label nCmps,
                void* databuf,
                char* datasetName,
                ioScalar atype
            );


        // Functions for handling mesh

            //- Mesh define wrapper
            void meshDefine(label regionID);

            //- Define points
            void meshDefinePoints(const fvMesh& m, label regionID);

            //- Define cells
            void meshDefineCells(const fvMesh& m, label regionID);

            //- Define Boundary meshes
            void meshDefineBoundaries(const fvMesh& m, label regionID);

            //- Mesh write wrapper
            void meshWrite(label regionID);

            //- Write points
            void meshWritePoints(const fvMesh& m, label regionID);

            //- Write cells
            void meshWriteCells(const fvMesh& m, label regionID);


        // Miscallenous definitions

            //- Read info for one region from it's own dictionary
            void read_region(const dictionary& dict, regionInfo& r);

            //- Append fieldName to the appropriate group
            label appendFieldGroup(regionInfo& r, const word& fieldName, const word& fieldType);

            //- Classify field types, returns the number of fields
            label classifyFields();
            //void test_print_obr();

            //- Disallow default bitwise copy construct
            adiosWrite(const adiosWrite&);

            //- Disallow default bitwise assignment
            void operator=(const adiosWrite&);

            //- Define ADIOS variables (before writing them)
            void defineVars();
            void deleteDefinitions();

        //- Pointers to cell shape models
        static const cellModel* unknownModel;
        static const cellModel* tetModel;
        static const cellModel* pyrModel;
        static const cellModel* prismModel;
        static const cellModel* hexModel;
        static const cellModel* wedgeModel;


        // Restart functions

            bool readData(const fileName&) const;

            bool readData(const instant& when) const;

            //- Read data from a given time
            // attempt to read data for all variables existing in memory
            // returns true on success, false on failure
            bool readData() const;
            bool readScalarFields(label regionID) const;
            bool readClouds(label regionID) const;
            bool readADIOSVar
            (
                ADIOS_FILE *f, ADIOS_SELECTION *sel,
                const char *path1, const char *path2,
                void *data
            ) const;


public:

    //- Runtime type information
    TypeName("adiosWrite");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        adiosWrite
        (
            const word& groupName,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~adiosWrite();


    // Member Functions

        //- Read the adiosWrite data
        virtual void read(const dictionary&);

        //- Execute, handles restart at first call if requested
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++,
        //  Note: currently does nothing
        virtual void timeSet();

        //- Write a timestep to file
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "adiosWriteFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
