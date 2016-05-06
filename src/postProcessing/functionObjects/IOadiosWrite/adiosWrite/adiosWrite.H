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

/**
 * write with cell-shapes to XDMF ids - will fail for non-primitive cell shapes
 */
// #define FOAM_ADIOS_CELL_SHAPES
// #undef  FOAM_ADIOS_CELL_SHAPES

/**
 * write internal fields and patch fields separately
 * poor prospects for restart
 */
// #define FOAM_ADIOS_PATCH_WRITE
// #undef  FOAM_ADIOS_PATCH_WRITE

/**
 * write cloud fields via explicit definition
 * poor prospects for restart
 */
// #define FOAM_ADIOS_CLOUD_EXPLICIT_NAMES
// #undef  FOAM_ADIOS_CLOUD_EXPLICIT_NAMES

/**
 * expand the particle information into separate fields
 * useful for debugging
 */
// #define FOAM_ADIOS_CLOUD_EXPAND
// #undef  FOAM_ADIOS_CLOUD_EXPAND

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// If OpenFOAM is compiled with single precision, the IO will be single
// precision.
#if defined(WM_SP) || defined(WRITE_SP)
# define ADIOS_IOSCALAR adios_real
typedef float ioScalar;
#elif defined(WM_DP)
# define ADIOS_IOSCALAR adios_double
typedef double ioScalar;
#else
# error "Neither WM_DP nor WM_SP defined when compiling adiosCore"
#endif

// Forward declaration of classes
class objectRegistry;
class dictionary;

/*---------------------------------------------------------------------------*\
                         Class adiosWrite Declaration
\*---------------------------------------------------------------------------*/

class adiosWrite
:
    public adiosCore
{
private:
    // Demand driven private data
    mutable Map<label>* shapeLookupPtr_;

    //- Is field-type supported?
    static bool supportedFieldType(const word& fieldType);


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
            //- Normally want auto-write?
            Switch autoWrite_;

            //- Name of the region
            word name_;

            //- Index of the region in the list
            label index_;

            //- Names of volume fields per region to write
            wordReList objectNames_;

            //- Names of clouds to write
            wordReList cloudNames_;

            //- Which (additional) attributes to write for particle data
            wordReList cloudAttribs_;

            //- Categorized scalar/vector/tensor vol fields
            fieldGroup<scalar> scalarFields_;
            fieldGroup<vector> vectorFields_;

            //- Categorized scalar/vector/tensor surf fields
            fieldGroup<scalar> surfaceScalarFields_;

#ifdef FOAM_ADIOS_CELL_SHAPES
            //- Number of data points in output Cell array on each process
            //  calculated in meshDefineCells, used in meshWriteCells
            //  only the local size is used currently
            List<label> cellDataSizes_;
#endif

            //- Array containing number of particles per process (often used list)
            //  Note this will have issues when multiple clouds exist
            List<label> nParticles_;

            //- Total number of particles (sum of nParticles_ list)
            //  Note this will have issues when multiple clouds exist
            label nTotalParticles_;

            //- Construct null
            regionInfo()
            :
                autoWrite_(true)
            {};

            ~regionInfo() {};


            //- If auto-write is in effect
            bool autoWrite() const
            {
                return autoWrite_;
            }


            //- Update the auto-write value if found in the dictionary
            bool getAutoWrite(const dictionary& dict)
            {
                return autoWrite_.readIfPresent("autoWrite", dict);
            }


            //- Read info for one region from its own dictionary
            void read
            (
                const fvMesh& mesh,
                const dictionary& dict,
                const dictionary& topDict = dictionary::null
            );

            //- Clear field groups
            void clearFields()
            {
                scalarFields_.clear(); // clear it because we will add all of them again and again
                vectorFields_.clear();
                surfaceScalarFields_.clear();
            }

            //- Append fieldName to the appropriate group
            label appendFieldGroup
            (
                const word& fieldName,
                const word& fieldType
            );

            //- Classify fields into the appropriate groups,
            //  returns the number of fields
            label classifyFields(const fvMesh& mesh);


            //- Path name for regions
            fileName regionPath() const
            {
                return "region" + Foam::name(index_);
            }

            //- Path name for fields
            fileName fieldPath(const string& var = string::null) const
            {
                return regionPath() / "field" / var;
            }

            //- Path name for meshes
            fileName meshPath(const string& var = string::null) const
            {
                return regionPath() / "polyMesh" / var;
            }

            //- Path name for cloud-Id
            fileName cloudPath(label cloudId) const
            {
                return regionPath() / ("cloud" + Foam::name(cloudId));
            }

            //- Simple text representation for output purposes
            word info() const
            {
                return "region" + Foam::name(index_) + "=" + name_;
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

        //- ADIOS input read method (default: BP)
        word readMethod_;

        //- ADIOS output write method (default: MPI)
        word writeMethod_;

        //- ADIOS output write parameters (default: "")
        string writeParams_;

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

protected:

    // Protected Member Functions

        // General functions

            //- Create and open dataset for all data
            void open();

            //- Close dataset
            void close();

            //- Mapping from OpenFOAM cell shapes to XDMF
            const Map<label>& shapeLookupMap() const;


        // General adios handling

            //- Define a variable with count elements
            //  Increment the output size and return the storage size
            size_t defineVariable
            (
                const char* name,
                enum ADIOS_DATATYPES type,
                size_t count
            );

            //- Define a variable with count elements
            //  Increment the output size and return the storage size
            size_t defineVariable
            (
                const string& name,
                enum ADIOS_DATATYPES type,
                size_t count
            );

            //- Define an integer variable containing a single element
            //  To be stored as part of a global 1D array of this name
            //  Increment the output size and return the storage size
            size_t defineIntVariable
            (
                const char* name
            );

            //- Define an integer variable containing a single element
            //  To be stored as part of a global 1D array of this name
            //  Increment the output size and return the storage size
            size_t defineIntVariable
            (
                const string& name
            );

            //- Define an integer variable with count elements
            //  Increment the output size and return the storage size
            size_t defineIntVariable
            (
                const char* name,
                size_t count
            );

            //- Define an integer variable with count elements
            //  Increment the output size and return the storage size
            size_t defineIntVariable
            (
                const string& name,
                size_t count
            );

            //- Define a double variable containing a single element
            //  To be stored as part of a global 1D array of this name
            //  Increment the output size and return the storage size
            size_t defineScalarVariable
            (
                const char* name
            );

            //- Define a double variable containing a single element
            //  To be stored as part of a global 1D array of this name
            //  Increment the output size and return the storage size
            size_t defineScalarVariable
            (
                const string& name
            );

            //- Define a double variable with count elements
            //  Increment the output size and return the storage size
            size_t defineScalarVariable
            (
                const char* name,
                size_t count
            );

            //- Define a double variable with count elements
            //  Increment the output size and return the storage size
            size_t defineScalarVariable
            (
                const string& name,
                size_t count
            );

            //- Define a byte-stream variable with count elements
            //  Increment the output size and return the storage size
            size_t defineStreamVariable
            (
                const char* name,
                size_t count
            );

            //- Define a byte-stream variable with count elements
            //  Increment the output size and return the storage size
            size_t defineStreamVariable
            (
                const string& name,
                size_t count
            );

            //- Define a vector (3 items) variable containing a single element
            //  Increment the output size and return the storage size
            size_t defineVectorVariable
            (
                const char* name
            );

            //- Define a vector (3 items) variable containing a single element
            //  Increment the output size and return the storage size
            size_t defineVectorVariable
            (
                const string& name
            );

            //- Define a vector (3 items) variable with count elements
            //  Increment the output size and return the storage size
            size_t defineVectorVariable
            (
                const char* name,
                size_t count
            );

            //- Define a vector (3 items) variable with count elements
            //  Increment the output size and return the storage size
            size_t defineVectorVariable
            (
                const string& name,
                size_t count
            );

            //- Define a string attribute
            void defineAttribute
            (
                const char* attrName,
                const char* varName,
                const string& value
            );

            //- Define a string attribute
            void defineAttribute
            (
                const char* attrName,
                const string& varName,
                const string& value
            );

            //- Define an int attribute
            void defineIntAttribute
            (
                const char* attrName,
                const char* varName,
                const int value
            );

            //- Define an int attribute
            void defineIntAttribute
            (
                const char* attrName,
                const string& varName,
                const int value
            );

            //- Define an integer array attribute
            bool defineListAttribute
            (
                const char* attrName,
                const string& varName,
                const UList<int>& list
            );

            //- Define a string array attribute
            bool defineListAttribute
            (
                const char* attrName,
                const string& varName,
                const wordList& list
            );

            //- Define a string array attribute
            bool defineListAttribute
            (
                const char* attrName,
                const string& varName,
                const stringList& list
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
                const string& name,
                const void* value
            );

            //- Write a variable from the list contents
            template<class T>
            void writeVariable
            (
                const char* name,
                const UList<T>& list
            )
            {
                writeVariable(name, list.cdata());
            }

            //- Write a variable from the list contents
            template<class T>
            void writeVariable
            (
                const string& name,
                const UList<T>& list
            )
            {
                writeVariable(name, list.cdata());
            }

            //- Write a variable with a single integer
            void writeIntVariable
            (
                const char* name,
                const int value
            );

            //- Write a variable with a single integer
            void writeIntVariable
            (
                const string& name,
                const int value
            );

            //- Write a variable with a single double value
            void writeScalarVariable
            (
                const char* name,
                const double value
            );

            //- Write a variable with a single double value
            void writeScalarVariable
            (
                const string& name,
                const double value
            );

            //- Write a variable with a single vector value
            void writeVectorVariable
            (
                const char* name,
                const vector& value
            );

            //- Write a variable with a single vector value
            void writeVectorVariable
            (
                const string& name,
                const vector& value
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
            size_t cloudDefine(regionInfo&);

            //- Cloud write wrapper
            void cloudWrite(const regionInfo&);

            //- Cloud define wrapper
            size_t cloudDefineExplicit(regionInfo&);

            //- Cloud write wrapper
            void cloudWriteExplicit(const regionInfo&);

            //- Write attributes for cloud particles
            void cloudWriteAttrib
            (
                label myParticles,
                label offset,
                label nTot,
                label nCmps,
                void* databuf,
                char* datasetName,
                ioScalar attrType
            );


        // Functions for handling mesh

            //- Mesh define wrapper
            size_t meshDefine(regionInfo&);

            //- Mesh write wrapper
            void meshWrite(const regionInfo&);


#ifdef FOAM_ADIOS_CELL_SHAPES

            //- Define variable and size for cell shapes
            void meshDefineCellShapes(const fvMesh&, regionInfo&);

            //- Write cell shapes
            void meshWriteCellShapes(const fvMesh&, const regionInfo&);
#endif


        // Miscallenous definitions

            //- Classify field types, returns the number of fields
            label classifyFields();

            //- Disallow default bitwise copy construct
            adiosWrite(const adiosWrite&);

            //- Disallow default bitwise assignment
            void operator=(const adiosWrite&);

            //- Define ADIOS variables (before writing them)
            //  Optionally define mesh information too
            size_t defineVars(const bool updateMesh);

            void deleteDefinitions();

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
