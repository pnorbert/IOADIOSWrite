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

Class
    Foam::adiosReader

Description

SourceFiles
    adiosReader.C

\*---------------------------------------------------------------------------*/

#ifndef adiosReader_H
#define adiosReader_H

#include "fileName.H"
#include "IBufStream.H"
#include "DynamicCharList.H"
#include "HashTable.H"

#include "Ostream.H"

#include "adios.h"
#include "adios_read.h"

#include "adiosTime.H"
#include "ParticleBinaryBlob.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class adiosReader Declaration
\*---------------------------------------------------------------------------*/

class adiosReader
{
public:

        // Forward declarations of various sub-components
        class VarInfo;
        class fieldInfo;
        class cloudInfo;

        //- Container type for attributes
        typedef HashTable<int, fileName> AttrContainer;

        //- Container type for variables
        typedef HashTable<VarInfo, fileName> VarContainer;

private:

    // Private Data

        ADIOS_FILE *file;
        ADIOS_SELECTION *selection;
        size_t maxLen;

        // simple hash of attributes names and their internal position (on read)
        AttrContainer attributes;

        VarContainer variables;

        wordList regionNames_;

        HashTable<wordList> cloudNames_;   // per region


    // Private Member Functions

        //- Disallow default bitwise copy construct
        adiosReader(const adiosReader&);

        //- Disallow default bitwise assignment
        void operator=(const adiosReader&);

        //- scan opened file for variables/attributes
        void scan(bool verbose=false);

        void select(ADIOS_SELECTION *sel);

public:

    // Constructors

        //- Null constructor
        adiosReader();

        //- Open specified filename and scan variables
        adiosReader(const fileName& bpFile, MPI_Comm comm);

    //- Destructor
        ~adiosReader();


        //- Open specified filename, select local process and scan variables
        bool open(const fileName& bpFile, MPI_Comm comm);

        //- Close file and reset
        void close();

        //- True if the file has been opened
        bool isGood() const;


        //- The max buffer size for the process local variables
        inline size_t sizeOf() const;


        HashTable<fieldInfo> getFieldInfo(const word& regName) const;
        HashTable<cloudInfo> getCloudInfo(const word& regName) const;


        //- True if the specified attribute exists
        inline bool hasAttribute(const string& varName) const;

        //- True if the specified variable exists
        inline bool hasVariable(const string& varName) const;


        //- Read integer attribute if it exists
        bool readIntAttributeIfPresent
        (
            const string& attrName,
            label& value
        ) const;


        //- Read string-list attribute if it exists
        bool readIntListAttributeIfPresent
        (
            const string& attrName,
            List<label>& lst
        ) const;


        //- Read string attribute if it exists
        bool readStringAttributeIfPresent
        (
            const string& attrName,
            string& value
        ) const;


        //- Read string-list attribute if it exists
        template<class StringType>
        bool readStringListAttributeIfPresent
        (
            const string& attrName,
            List<StringType>& lst
        ) const;


        //- Read integer variable (single-value) if it exists
        bool readIntVariableIfPresent
        (
            const string& varName,
            label& value
        ) const;


        //- Read scalar variable (single-value) if it exists
        bool readScalarVariableIfPresent
        (
            const string& varName,
            scalar& value
        ) const;


        //- Get a mandatory integer attribute
        inline label getIntAttribute(const string& name) const;

        //- Get a mandatory integer list attribute
        inline List<label> getIntListAttribute(const string& name) const;

        //- Get a mandatory string attribute
        inline string getStringAttribute(const string& name) const;

        //- Get a mandatory integer attribute
        template<class StringType>
        inline List<StringType> getStringListAttribute(const string& name) const;

        //- Get a mandatory integer variable (single-value)
        inline label getIntVariable(const string& name) const;

        //- Get a mandatory scalar variable (single-value)
        inline scalar getScalarVariable(const string& name) const;


        //- Retrieve the variable directly into the memory specified
        bool getVariable(const string& varName, void* data) const;


        //- Retrieve the variable into our buffer area
        //  On success the buffer is resized according to the data.
        //  On failure the buffer is resized to zero.
        template<class T>
        size_t getBuffered(const string& varName, DynamicList<T>& buffer) const;

};


//- Partial transcription of ADIOS variable information (ADIOS_VARINFO)
class adiosReader::VarInfo
{
    // Private data

        fileName name_;
        label varid_;
        enum ADIOS_DATATYPES type_;

        label nElem_;
        size_t nBytes_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const VarInfo&);

        //- Read variable information from bp file
        bool readInfo(ADIOS_FILE *file, bool verbose=false);

public:

    // Constructors

        //- Construct empty
        VarInfo(const char* varName);

        //- Construct empty
        VarInfo(const string& varName);

        //- Construct by reading variable information from bp file
        VarInfo(ADIOS_FILE *file, const char* varName, bool verbose=false);

        //- Construct by reading variable information from bp file
        VarInfo(ADIOS_FILE *file, const string& varName, bool verbose=false);


    //- Destructor
    ~VarInfo();


    // Member Functions

        //- The variable name
        inline const fileName& name() const;

        //- The underlying ADIOS data type
        inline enum ADIOS_DATATYPES dataType() const;

        //- The process local number of elements (as flattened dimensions)
        inline label nElem() const;

        //- The byte-size for the process local elements
        inline size_t sizeOf() const;

};


class adiosReader::fieldInfo
:
    public fileName
{
    word className_;
    size_t nBytes_;

public:

    fieldInfo()
    :
        fileName(),
        className_(),
        nBytes_(0)
    {}

    fieldInfo
    (
        const fileName& varName,
        size_t nbytes,
        const word& className
    )
    :
        fileName(varName),
        className_(className),
        nBytes_(nbytes)
    {}


    //- The region name is before the first /
    const word regionName() const
    {
        return substr(0, find('/')); // not designed for failure
    }


    //- Return the full variable name
    inline const fileName& fullName() const
    {
        return static_cast<const fileName&>(*this);
    }


    //- Return file name (part beyond last /)
    using fileName::name;

    //- Return directory path name (part before last /)
    using fileName::path;


    //- The class (volScalarField etc)
    inline const word& type() const
    {
        return className_;
    }

    //- Local size in bytes
    inline size_t sizeOf() const
    {
       return nBytes_;
    }
};


// currently looks mostly identical to fieldInfo, but this will change
class adiosReader::cloudInfo
:
    public fileName
{
    word className_;
    size_t nBytes_;

    //- Total number of particles (sum of nParticle_ list)
    label nTotal_;

    //- Number of particles on local process
    label nParticle_;

    //- blob width
    label width_;

    //- blob width
    ParticleBinaryBlob blob_;


    bool read(const adiosReader&);


public:


    cloudInfo();
    cloudInfo(const fileName& varName);
    cloudInfo
    (
        const adiosReader::VarInfo& varinfo,
        const adiosReader& reader
    );


    //- The region name is before the first /
    const word regionName() const
    {
        return substr(0, find('/')); // not designed for failure
    }


    //- Return the full variable name
    inline const fileName& fullName() const
    {
        return static_cast<const fileName&>(*this);
    }


    //- Return file name (part beyond last /)
    using fileName::name;

    //- Return directory path name (part before last /)
    using fileName::path;


    //- The class (volScalarField etc)
    inline const word& type() const
    {
        return className_;
    }

    //- Total number of particles across all processes
    label nTotal() const
    {
        return nTotal_;
    }


    //- The number of particles on this process
    label nParticle() const
    {
        return nParticle_;
    }


    //- Local size in bytes
    inline size_t sizeOf() const
    {
        return nBytes_;
    }


    //- Blob width in bytes
    inline size_t width() const
    {
        return width_;
    }


    //- The blob information
    const ParticleBinaryBlob& blob() const
    {
        return blob_;
    }


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "adiosReaderI.H"

#ifdef NoRepository
    #include "adiosReaderTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
