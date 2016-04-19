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
#include "adiosReader.H"

#include "runTimeSelectionTables.H"
#include "interpolation.H"
#include "fvMesh.H"
#include "surfaceMesh.H"
#include "fvsPatchFields.H"

#include "adios.h"
#include "adios_read.h"
#include <stdio.h> // sprintf

#define FOAM_ADIOS_PATCH_WRITE
// #undef  FOAM_ADIOS_PATCH_WRITE

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// If OpenFOAM is compiled with single precision, the IO will be single
// precision.
#if defined(WM_SP) || defined(WRITE_SP)
# define ADIOS_SCALAR   adios_real
typedef float ioScalar;
#elif defined(WM_DP)
# define ADIOS_SCALAR   adios_double
typedef double ioScalar;
#else
# error "Neither WM_DP nor WM_SP defined when compiling adiosCore"
#endif

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
            regionInfo() {};
            ~regionInfo() {};

            //- Name of the region
            word name_;

            //- Index of the region in the list
            label index_;

            //- Names of volume fields per region to write
            wordReList objectNames_;

            //- Names of clouds to write
            wordReList cloudNames_;

            //- Which attributes to write for particle data
            wordReList cloudAttribs_;

            //- Categorized scalar/vector/tensor vol fields
            fieldGroup<scalar> scalarFields_;
            fieldGroup<vector> vectorFields_;

            //- Categorized scalar/vector/tensor surf fields
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

            //- The standard path name for a field
            fileName fieldVarPath(const std::string& var) const
            {
                return fileName
                (
                    "region" + Foam::name(index_) / "field" / var
                );
            }

            //- The standard path name for a mesh variable
            fileName meshVarPath(const std::string& var) const
            {
                return fileName
                (
                    "mesh" + Foam::name(index_) / var
                );
            }
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

        //- Buffer for reading/writing streamed content to/from adious files
        DynamicCharList iobuffer_;

    // Protected Member Functions

        // General functions

            // Create and open dataset for all data
            void open();

            // Close dataset
            void close();


        // General adios handling

            //- Define a variable with count elements
            //  Increment the output size and return the variable size
            size_t defineVariable
            (
                const char* name,
                enum ADIOS_DATATYPES type,
                size_t count = 1
            );

            //- Define a variable with count elements
            //  Increment the output size and return the variable size
            size_t defineVariable
            (
                const fileName& name,
                enum ADIOS_DATATYPES type,
                size_t count = 1
            );

            //- Define a vector (3 items) variable with count elements
            //  Increment the output size and return the variable size
            size_t defineVectorVariable
            (
                const char* name,
                enum ADIOS_DATATYPES type,
                size_t count
            );

            //- Define a vector (3 items) variable with count elements
            //  Increment the output size and return the variable size
            size_t defineVectorVariable
            (
                const fileName& name,
                enum ADIOS_DATATYPES type,
                size_t count
            );

            //- Define a string attribute
            void defineAttribute
            (
                const char* attrName,
                const char* varName,
                const std::string& value
            );

            //- Define a string attribute
            void defineAttribute
            (
                const char* attrName,
                const fileName& varName,
                const std::string& value
            );

            //- Define a bool attribute (save as byte)
            void defineAttribute
            (
                const char* attrName,
                const char* varName,
                const bool value
            );

            //- Define a bool attribute (save as byte)
            void defineAttribute
            (
                const char* attrName,
                const fileName& varName,
                const bool value
            );

            //- Write a variable
            void writeVariable
            (
                const char* name,
                const void* value
            );

            //- Write a variable
            void writeVariable
            (
                const fileName& name,
                const void* value
            );


        // Functions for handling field data

            //- Field define wrapper
            //  Return the maximum associated stream length
            size_t fieldDefine(const regionInfo&);

            //- Field define for a given type
            //  Return the maximum associated stream length
            template<class FieldType>
            size_t fieldDefine
            (
                const fvMesh& mesh,
                const regionInfo&,
                const fieldGroup<typename FieldType::value_type>& fields
            );

            //- Field write wrapper
            void fieldWrite(const regionInfo&);

            //- Write data for a given type
            template<class FieldType>
            void fieldWrite
            (
                const fvMesh& mesh,
                const regionInfo&,
                const fieldGroup<typename FieldType::value_type>& fields
            );

            //- Field read wrapper
            bool fieldRead
            (
                adiosReader::helper& helper,
                regionInfo&
            );

            //- Read data for a given type
            template<class FieldType>
            bool fieldRead
            (
                adiosReader::helper& helper,
                const fvMesh& mesh,
                regionInfo&,
                const fieldGroup<typename FieldType::value_type>& fields
            );


      // Functions for handling particle cloud data

            //- Cloud define wrapper
            void cloudDefine(regionInfo&);

            //- Cloud write wrapper
            void cloudWrite(const regionInfo&);

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
            size_t meshDefine(regionInfo&);

            //- Define points
            size_t meshDefinePoints(const fvMesh&, regionInfo&);

            //- Define cells
            size_t meshDefineCells(const fvMesh&, regionInfo&);

            //- Define Boundary meshes
            size_t meshDefineBoundaries(const fvMesh&, regionInfo&);

            //- Mesh write wrapper
            void meshWrite(const regionInfo&);

            //- Write points
            void meshWritePoints(const fvMesh&, const regionInfo&);

            //- Write cells
            void meshWriteCells(const fvMesh&, const regionInfo&);


        // Miscallenous definitions

            //- Read info for one region from it's own dictionary
            void read_region(const dictionary& dict, regionInfo&);

            //- Append fieldName to the appropriate group
            label appendFieldGroup(regionInfo&, const word& fieldName, const word& fieldType);

            //- Classify field types, returns the number of fields
            label classifyFields();
            //void test_print_obr();

            //- Disallow default bitwise copy construct
            adiosWrite(const adiosWrite&);

            //- Disallow default bitwise assignment
            void operator=(const adiosWrite&);

            //- Define ADIOS variables (before writing them)
            //  Optionally define mesh information too
            size_t defineVars(const bool updateMesh);

            void deleteDefinitions();

        //- Pointers to cell shape models
        static const cellModel* unknownModel;
        static const cellModel* tetModel;
        static const cellModel* pyrModel;
        static const cellModel* prismModel;
        static const cellModel* hexModel;
        static const cellModel* wedgeModel;


        // Restart functions

            bool readData(const fileName&);
            bool readData(const instant& when);

            //- Read data from a given time
            // attempt to read data for all variables existing in memory
            // returns true on success, false on failure
            bool readData();

            bool readClouds(adiosReader::helper&, regionInfo&);

            // bool readScalarFields(adiosReader::helper&, regionInfo&);
            // NOT YET: bool readClouds(adiosReader::helper&, label regionID);

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
    #include "adiosReadFieldTemplates.C"
    #include "adiosWriteFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
