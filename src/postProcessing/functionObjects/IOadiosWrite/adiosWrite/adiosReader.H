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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class adiosReader Declaration
\*---------------------------------------------------------------------------*/

class adiosReader
{
public:

    class helper
    {

        // Private Member Functions

            //- Disallow default bitwise copy construct
            helper(const helper&);

            //- Disallow default bitwise assignment
            void operator=(const helper&);

    public:
        DynamicCharList& buffer;
        ADIOS_FILE *file;
        ADIOS_SELECTION *selection;
        size_t maxLen;

        // simple hash of attributes names and their internal position (on read)
        HashTable<int, fileName> attributes;

        HashTable<size_t, fileName> variables;

        wordList regionNames_;
        HashTable<wordList> cloudNames_;   // per region

        helper(const DynamicCharList&);

        ~helper();


        bool open(const fileName& bpFile, MPI_Comm comm);

        void select(ADIOS_SELECTION *sel);

        //- Return the number of bytes associated with the dataset
        //  Values are for the current processor.
        size_t sizeOf(const char* datasetName, bool verbose=false);

        //- Return the number of bytes associated with the dataset
        //  Values are for the current processor.
        size_t sizeOf(const string& datasetName, bool verbose=false);

        //- Retrieve named dataset into our buffer area
        //  On success the buffer is resized according to the data.
        //  On failure the buffer is resized to zero.
        bool getDataSet(const string& datasetName);

        //- Retrieve named dataset directly into the memory specified
        bool getDataSet(const string& datasetName, void* data);

        void scan(bool verbose=false);

        void close();

        //- Read integer attribute if it exists
        bool readIntAttributeIfPresent
        (
            const string& attrName,
            label &value
        ) const;


        //- Read string-list attribute if it exists
        bool readIntListAttributeIfPresent
        (
            const string& attrName,
            List<label>& lst
        ) const;


        //- Get a mandatory integer attribute
        label getIntAttribute(const string& attrName) const;


        //- Get a mandatory integer list attribute
        List<label> getIntListAttribute(const string& attrName) const;


        //- Read string-list attribute if it exists
        template<class StringType>
        bool readStringListAttributeIfPresent
        (
            const string& attrName,
            List<StringType>& lst
        ) const;


        //- Get a mandatory integer attribute
        template<class StringType>
        List<StringType> getStringListAttribute(const string& attrName) const;

    };

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "adiosReaderTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
