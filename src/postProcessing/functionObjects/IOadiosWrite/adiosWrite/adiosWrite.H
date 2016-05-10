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
#include "adiosTime.H"

#include "runTimeSelectionTables.H"
#include "interpolation.H"
#include "fvMesh.H"
#include "surfaceMesh.H"
#include "fvsPatchFields.H"

#include "adios.h"
#include "adios_read.h"

/**
 * expand the particle information into separate fields
 * useful for debugging
 */
// #define FOAM_ADIOS_CLOUD_EXPAND
// #undef  FOAM_ADIOS_CLOUD_EXPAND

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

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

    //- Is field-type supported?
    static bool supportedFieldType(const word& fieldType);


protected:

    // Protected classes

        // Per-region info variables grouped together to be able to create a list of them
        class cloudInfo
        {
            //- Name of the cloud
            word name_;

            //- Total number of particles (sum of nParticle_ list)
            label nTotal_;

            //- Array containing number of particles per process (often used list)
            List<label> nParticle_;

        public:

            //- Construct null
            cloudInfo(const word& name)
            :
               name_(name),
               nTotal_(0),
               nParticle_(Pstream::nProcs())
            {
                reset();
            }

            ~cloudInfo() {}


            //- Name of the cloud
            const word& name() const
            {
                return name_;
            }


            //- Total number of particles across all processes
            label nTotal() const
            {
                return nTotal_;
            }


            //- The number of particles on this process
            label nParticle() const
            {
                return nParticle_[Pstream::myProcNo()];
            }


            //- Set the number of particles and sync across processes
            //  Return the overall total
            label nParticle(label n)
            {
                // The number of particles on each process
                nParticle_[Pstream::myProcNo()] = n;
                Pstream::gatherList(nParticle_);
                Pstream::scatterList(nParticle_);

                // Sum total number of particles on all processes
                nTotal_ = sum(nParticle_);

                return nTotal_;
            }


            //- Local process offset within a global array of parcels
            label offset() const
            {
                label off = 0;
                for (label prev=0; prev < Pstream::myProcNo(); ++prev)
                {
                    off += nParticle_[prev];
                }

                return off;
            }


            //- Reset to zeros
            void reset()
            {
                nTotal_ = 0;
                forAll(nParticle_, i)
                {
                    nParticle_[i] = 0;
                }
            }

        };


        // Per-region info variables grouped together to be able to create a list of them
        class regionInfo
        {
        public:
            //- Normally want auto-write?
            Switch autoWrite_;

            //- Name of the region
            word name_;

            //- Names of volume fields per region to write
            wordReList objectNames_;

            //- Names of clouds to write
            wordReList cloudNames_;

            //- Which (additional) attributes to write for particle data
            wordReList cloudAttribs_;

            //- Fields to be written. Field name as key, class as value
            HashTable<word> fieldsToWrite_;

            //- Clouds to be processed
            SLList<cloudInfo> cloudInfo_;

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

            //- Clear field to be written
            void clearFields()
            {
                fieldsToWrite_.clear();
            }

            //- Classify fields into the appropriate groups,
            //  returns the number of fields
            label classifyFields
            (
                const fvMesh& mesh,
                bool verbose = false
            );


            //- Path name for regions
            string regionPath() const
            {
                return adiosCore::regionPath(name_);
            }

            //- Path name for fields
            fileName fieldPath(const string& var) const
            {
                return adiosCore::fieldPath(name_, var);
            }


            //- Path name for named cloud
            fileName cloudPath(const string& cloudName) const
            {
                return adiosCore::cloudPath(name_, cloudName);
            }


            //- Path name for named cloud
            fileName cloudPath(const cloudInfo& cinfo) const
            {
                return cloudPath(cinfo.name());
            }


            //- Simple text representation for output purposes
            word info() const
            {
                return name_;
            }

        };


    // Protected data

        // Reference to the object registry
        const objectRegistry& obr_;

        //- Reference to the primary mesh database
        const fvMesh& primaryMesh_; //  primaryMesh_(refCast<const fvMesh>(obr))

        //- Reference to the time database
        const Time& time_;  // = primaryMesh_.time()

        //- track time index from restart (avoid write immediately after restart)
        int restartIndex_;

        // Restart time requested
        scalar restartTime_;

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
                const FieldType& field,
                const fileName& varPath
            );


            //- Field define for a given type
            //  Return the maximum associated stream length
            template<class FieldType>
            size_t fieldDefineInternal
            (
                const FieldType& field,
                const fileName& varPath
            );

            //- Field write wrapper
            void fieldWrite(const regionInfo&);

            //- Write data for a given type
            template<class FieldType>
            void fieldWrite
            (
                const FieldType& field,
                const fileName& varPath
            );

            //- Read data for a given volume field type
            bool readVolField
            (
                regIOobject* obj,
                const adiosReader& reader,
                const adiosReader::fieldInfo& src
            );

            //- Read data for a given volume field type
            bool readCloud
            (
                const fvMesh& mesh,
                const word& cloudName,
                const adiosReader& reader,
                const adiosReader::cloudInfo& src
            );

            //- Read data for a given type
            template<class FieldType>
            bool fieldRead
            (
                FieldType& field,
                const adiosReader& reader,
                const adiosReader::fieldInfo& src
            );

            //- Read data for a given type
            template<class FieldType>
            bool fieldRead
            (
                const adiosReader& reader,
                const fvMesh& mesh,
                const adiosReader::fieldInfo& src
            );


      // Functions for handling particle cloud data

            //- Define ADIOS cloud variables/attributes
            size_t cloudDefine(regionInfo&);

            //- Write ADIOS cloud variables
            void cloudWrite(const regionInfo&);

        // Functions for handling mesh

            //- Define ADIOS mesh variables
            size_t meshDefine
            (
                const fvMesh&,
                const bool updateMesh
            );

            //- Write ADIOS mesh variables
            void meshWrite
            (
                const fvMesh&,
                const bool updateMesh
            );


        // Miscallenous definitions

            //- Classify field types, returns the number of fields
            label classifyFields
            (
                bool verbose = false
            );

            //- Disallow default bitwise copy construct
            adiosWrite(const adiosWrite&);

            //- Disallow default bitwise assignment
            void operator=(const adiosWrite&);

            //- Define ADIOS variables (before writing them)
            //  Optionally define mesh information too
            size_t defineVars(const bool updateMesh);

            void deleteDefinitions();

        // Restart functions

            adiosTime readData(const fileName&);
            adiosTime readData(const instant& when);

            //- Read data from a given time
            // attempt to read data for all variables existing in memory
            // returns true on success, false on failure
            adiosTime readData();

            //- Execute, handles restart at first call if requested
            bool restart();

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
    #include "adiosFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
