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

Class
    Foam::adiosCoreWrite

Description
    Base for writing ADIOS files.

SourceFiles
    adiosCoreWrite.C

\*---------------------------------------------------------------------------*/

#ifndef adiosCoreWrite_H
#define adiosCoreWrite_H

#include "adiosCore.H"
#include "dictionary.H"
#include "Pstream.H"
#include "DynamicCharList.H"

#include "adios.h"
#include "adios_read.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class TimeState;

/*---------------------------------------------------------------------------*\
                       Class adiosCoreWrite Declaration
\*---------------------------------------------------------------------------*/

class adiosCoreWrite
:
    public adiosCore
{
    // Private data

        //- ADIOS group identifier
        int64_t groupID_;

        //- ADIOS file handler (from open)
        int64_t fileID_;

        //- ADIOS input read method (default: BP)
        word readMethod_;

        //- ADIOS output write method (default: MPI)
        word writeMethod_;

        //- ADIOS output write parameters (default: "")
        string writeParams_;

        // MPI parallel code variables
        MPI_Comm comm_;  // Pstream does not provide communicator

        //- Map of adios variable names to their internal adios id
        typedef HashTable<int64_t, string> VarHashType;

        //- Map of adios variable names to their internal adios id
        VarHashType vars_;


    // Private Member Functions

        //- Read the adiosWriter data
        void read(const dictionary&);

        //- Disallow default bitwise copy construct
        adiosCoreWrite(const adiosCoreWrite&);

        //- Disallow default bitwise assignment
        void operator=(const adiosCoreWrite&);

protected:

    // Protected data

        //- Buffer for reading/writing streamed content to/from adios files
        DynamicCharList iobuffer_;


    // Static Member Functions

        //- Is field-type supported?
        static bool supportedFieldType(const word& fieldType);


    // Protected Member Functions

        //- MPI communicator
        MPI_Comm communicator() const
        {
            return comm_;
        }


    // General functions

        //- Create and open output file
        bool open(const fileName& dataFile);

        //- Close output file, flushing ADIOS buffer
        void close();

        //- Reset variable and attribute definitions
        void reset();

        //- The maximum size needed by any declared adios variable
        size_t maxSize() const;

        //- The output size for one timestep (Sum of declared variables).
        size_t outputSize() const;


    // Member Functions

        //- Define base attributes
        void defineBaseAttributes();

        //- Define attributes for OpenFOAM time information
        void defineTimeAttributes(const TimeState&);

        //- Define variables for OpenFOAM time information
        void defineTime();

        //- Write OpenFOAM time information
        void writeTime(const TimeState&);

        void definePatchAttributes(const polyMesh&);

        void defineMeshPoints(const polyMesh&);
        void defineMeshFaces(const polyMesh&);
        void defineMesh(const polyMesh&);
        void writeMeshPoints(const polyMesh&);
        void writeMeshFaces(const polyMesh&);
        void writeMesh(const polyMesh&);

        //- Define adios variable for field, return adios Id.
        template<class FieldType>
        int64_t defineField
        (
            const FieldType& field,
            const fileName& basePath
        );

        //- Define adios variable for internal field, return adios Id.
        template<class FieldType>
        int64_t defineInternalField
        (
            const FieldType& field,
            const fileName& basePath
        );

        //- Write field data for a given type
        template<class FieldType>
        void writeField
        (
            const FieldType& field,
            const fileName& basePath
        );


    // General adios handling

        //- Generic define a variable, return adios Id.
        int64_t defineVariable
        (
            const string& name,
            enum ADIOS_DATATYPES type,
            const string& localDimensions,
            const string& globalDimensions,
            const string& offsets
        );

        //- Define a variable with count elements, return adios Id.
        int64_t defineVariable
        (
            const char* name,
            enum ADIOS_DATATYPES type,
            size_t count
        );

        //- Define a variable with count elements, return adios Id.
        int64_t defineVariable
        (
            const string& name,
            enum ADIOS_DATATYPES type,
            size_t count
        );

        //- Define an integer variable containing a single element.
        //  To be stored as part of a global 1D array of this name.
        //  Return adios Id.
        int64_t defineIntVariable(const char* name);

        //- Define an integer variable containing a single element.
        //  To be stored as part of a global 1D array of this name.
        //  Return adios Id.
        int64_t defineIntVariable(const string& name);

        //- Define an integer variable with count elements.
        //  Return adios Id.
        int64_t defineIntVariable(const char* name, size_t count);

        //- Define an integer variable with count elements
        //  Return adios Id.
        int64_t defineIntVariable(const string& name, size_t count);

        //- Define a double variable containing a single element.
        //  To be stored as part of a global 1D array of this name.
        //  Return adios Id.
        int64_t defineScalarVariable(const char* name);

        //- Define a double variable containing a single element.
        //  To be stored as part of a global 1D array of this name.
        //  Return adios Id.
        int64_t defineScalarVariable(const string& name);

        //- Define a double variable with count elements.
        //  Return adios Id.
        int64_t defineScalarVariable(const char* name, size_t count);

        //- Define a double variable with count elements.
        //  Return adios Id.
        int64_t defineScalarVariable(const string& name, size_t count);

        //- Define a byte-stream variable with count elements.
        //  Return adios Id.
        int64_t defineStreamVariable(const char* name, size_t count);

        //- Define a byte-stream variable with count elements.
        //  Return adios Id.
        int64_t defineStreamVariable(const string& name, size_t count);

        //- Define a vector (3 items) variable containing a single element.
        //  Return adios Id.
        int64_t defineVectorVariable(const char* name);

        //- Define a vector (3 items) variable containing a single element.
        //  Return adios Id.
        int64_t defineVectorVariable(const string& name);

        //- Define a vector (3 items) variable with count elements.
        //  Return adios Id.
        int64_t defineVectorVariable(const char* name, size_t count);

        //- Define a vector (3 items) variable with count elements.
        //  Return adios Id.
        int64_t defineVectorVariable(const string& name, size_t count);


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

        //- Define a double attribute
        void defineScalarAttribute
        (
            const char* attrName,
            const char* varName,
            const double value
        );

        //- Define a double attribute
        void defineScalarAttribute
        (
            const char* attrName,
            const string& varName,
            const double value
        );

        //- Define an array attribute of integer values
        bool defineListAttribute
        (
            const char* attrName,
            const string& varName,
            const UList<int>& list
        );

        //- Define an array attribute of double values
        bool defineListAttribute
        (
            const char* attrName,
            const string& varName,
            const UList<double>& list
        );

        //- Define an array attribute of strings
        bool defineListAttribute
        (
            const char* attrName,
            const string& varName,
            const wordList& list
        );

        //- Define an array attribute of strings
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
        inline void writeVariable
        (
            const char* name,
            const UList<T>& list
        )
        {
            writeVariable(name, list.cdata());
        }

        //- Write a variable from the list contents
        template<class T>
        inline void writeVariable
        (
            const string& name,
            const UList<T>& list
        )
        {
            writeVariable(name, list.cdata());
        }

        //- Write a variable with 1 integer
        void writeIntVariable
        (
            const char* name,
            const int value
        );

        //- Write a variable with 1 integer
        void writeIntVariable
        (
            const string& name,
            const int value
        );

        //- Write a variable with 1 double value
        void writeScalarVariable
        (
            const char* name,
            const double value
        );

        //- Write a variable with 1 double value
        void writeScalarVariable
        (
            const string& name,
            const double value
        );

        //- Write a variable with 1 vector value
        void writeVectorVariable
        (
            const char* name,
            const vector& value
        );

        //- Write a variable with 1 vector value
        void writeVectorVariable
        (
            const string& name,
            const vector& value
        );

public:

    // Constructors

        //- Construct for specified group name, using dictionary parameters
        adiosCoreWrite(const word& groupName, const dictionary&);


    //- Destructor
    virtual ~adiosCoreWrite();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "adiosCoreWriteFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
