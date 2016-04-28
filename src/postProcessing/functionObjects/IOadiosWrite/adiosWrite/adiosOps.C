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

\*---------------------------------------------------------------------------*/

#include "adiosWrite.H"
#include "fileName.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::adiosWrite::open()
{
    Info<< "adiosWrite::open:" << endl;

    // Create output directory if initially non-existent
    static bool checkdir = true;
    if (checkdir && !isDir(dataDirectory))
    {
        mkDir(dataDirectory);
    }
    checkdir = false;


    // Find and create filename
    char mode[] = "w";

    fileName dataFile(dataDirectory/obr_.time().timeName() + ".bp");

    // Create/open a new file collectively.
    /*
    static int i = 0;
    if (timeSteps_ == 0)
    {
        mode[0] = 'w';
        do
        {
            //sprintf(dataFile, "adiosData/%s%i.bp", name_.c_str(), i);
            sprintf(dataFile, "%s%i.bp", name_.c_str(), i);
            i++;
        }
        while (isFile(dataFile));
    }
    else
    {
        mode[0] = 'a';
        sprintf(dataFile, "%s%i.bp", name_.c_str(), i);
    }
    */

    // Print info to terminal
    Info<< "  adiosWrite: Chosen filename " << dataFile << nl << endl;

    int err = adios_open
    (
        &fileID_,
        name().c_str(),     // group-name
        dataFile.c_str(),
        mode,
        comm_
    );

    if (!err)
    {
        // Tell ADIOS how many bytes we are going to write in this step (by this process)
        uint64_t total_size; // user data bytes + metadata, buffer size should be bigger
        adios_group_size(fileID_, outputSize_, &total_size);
    }
    else
    {
        WarningInFunction
            << "File " << name()
            << " could not be created: " << nl
            << adios_get_last_errmsg();

        fileID_ = 0;
    }
}


void Foam::adiosWrite::close()
{
    MPI_Barrier(MPI_COMM_WORLD);
    Info<< "adiosWrite::close" << nl << endl;

    // Close the file
    if (fileID_)
    {
        adios_close(fileID_);
        fileID_ = 0;
    }
}


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //


size_t Foam::adiosWrite::defineVariable
(
    const char* name,
    enum ADIOS_DATATYPES type,
    size_t count
)
{
    adios_define_var
    (
        groupID_,
        name,                           // name
        NULL,                           // path (deprecated)
        type,                           // data-type
        Foam::name(count).c_str(),      // local dimensions
        NULL,                           // global dimensions
        NULL                            // local offsets
    );

    int sz = count * adios_type_size(type, NULL);

    outputSize_ += sz;

    return sz;
}


size_t Foam::adiosWrite::defineVariable
(
    const fileName& name,
    enum ADIOS_DATATYPES type,
    size_t count
)
{
    return defineVariable(name.c_str(), type, count);
}


size_t Foam::adiosWrite::defineIntVariable
(
    const char* name
)
{
    adios_define_var
    (
        groupID_,
        name,                           // name
        NULL,                           // path (deprecated)
        adios_integer,                  // data-type
        "1",                            // local dimensions
        Foam::name(Pstream::nProcs()).c_str(),  // global 1D array of this info
        Foam::name(Pstream::myProcNo()).c_str() // offsets of this process into array
    );

    int sz = adios_type_size(adios_integer, NULL);

    outputSize_ += sz;

    return sz;
}


size_t Foam::adiosWrite::defineIntVariable
(
    const fileName& name
)
{
    return defineIntVariable(name.c_str());
}


size_t Foam::adiosWrite::defineIntVariable
(
    const char* name,
    size_t count
)
{
    adios_define_var
    (
        groupID_,
        name,                           // name
        NULL,                           // path (deprecated)
        adiosTraits<label>::adiosType,  // data-type
        Foam::name(count).c_str(),      // local dimensions
        NULL,                           // global dimensions
        NULL                            // local offsets
    );

    int sz = count * adiosTraits<label>::adiosSize;

    outputSize_ += sz;

    return sz;
}


size_t Foam::adiosWrite::defineIntVariable
(
    const fileName& name,
    size_t count
)
{
    return defineIntVariable(name.c_str(), count);
}


size_t Foam::adiosWrite::defineScalarVariable
(
    const char* name
)
{
    adios_define_var
    (
        groupID_,
        name,                           // name
        NULL,                           // path (deprecated)
        adiosTraits<scalar>::adiosType, // data-type
        "1",                            // local dimensions
        Foam::name(Pstream::nProcs()).c_str(),  // global 1D array of this info
        Foam::name(Pstream::myProcNo()).c_str() // offsets of this process into array
    );

    size_t sz = adiosTraits<scalar>::adiosSize;

    outputSize_ += sz;

    return sz;
}


size_t Foam::adiosWrite::defineScalarVariable
(
    const fileName& name
)
{
    return defineScalarVariable(name.c_str());
}


size_t Foam::adiosWrite::defineScalarVariable
(
    const char* name,
    size_t count
)
{
    adios_define_var
    (
        groupID_,
        name,                           // name
        NULL,                           // path (deprecated)
        adiosTraits<scalar>::adiosType, // data-type
        Foam::name(count).c_str(),      // local dimensions
        NULL,                           // global dimensions
        NULL                            // local offsets
    );

    int sz = count * adiosTraits<scalar>::adiosSize;

    outputSize_ += sz;

    return sz;
}


size_t Foam::adiosWrite::defineScalarVariable
(
    const fileName& name,
    size_t count
)
{
    return defineScalarVariable(name.c_str(), count);
}


size_t Foam::adiosWrite::defineStreamVariable
(
    const char* name,
    size_t count
)
{
    // use unsigned byte:
    // keeps people from thinking that the min/max statistics have any meaning

    adios_define_var
    (
        groupID_,
        name,                           // name
        NULL,                           // path (deprecated)
        adios_unsigned_byte,            // data-type
        Foam::name(count).c_str(),      // local dimensions
        NULL,                           // global dimensions
        NULL                            // local offsets
    );

    outputSize_ += count;

    return count;
}


size_t Foam::adiosWrite::defineStreamVariable
(
    const fileName& name,
    size_t count
)
{
    return defineStreamVariable(name.c_str(), count);
}


size_t Foam::adiosWrite::defineVectorVariable
(
    const char* name
)
{
    string globalDims = Foam::name(Pstream::nProcs()) + ",3";
    string offsetDims = Foam::name(Pstream::myProcNo());

    adios_define_var
    (
        groupID_,
        name,                           // name
        NULL,                           // path (deprecated)
        adiosTraits<scalar>::adiosType, // data-type
        "1,3",                          // local dimensions
        globalDims.c_str(),             // global dimensions
        offsetDims.c_str()              // offsets of this process into array
    );

    size_t sz = 3*adiosTraits<scalar>::adiosSize;

    outputSize_ += sz;

    return sz;
}


size_t Foam::adiosWrite::defineVectorVariable
(
    const fileName& name
)
{
    return defineVectorVariable(name.c_str());
}


size_t Foam::adiosWrite::defineVectorVariable
(
    const char* name,
    size_t count
)
{
    string localDims = Foam::name(count) + ",3";

    adios_define_var
    (
        groupID_,
        name,                           // name
        NULL,                           // path (deprecated)
        adiosTraits<scalar>::adiosType, // data-type
        localDims.c_str(),              // local dimensions
        NULL,                           // global dimensions
        NULL                            // local offsets
    );

    int sz = count * 3*adiosTraits<scalar>::adiosSize;
    outputSize_ += sz;

    return sz;
}


size_t Foam::adiosWrite::defineVectorVariable
(
    const fileName& name,
    size_t count
)
{
    return defineVectorVariable(name.c_str(), count);
}


void Foam::adiosWrite::defineAttribute
(
    const char* attrName,
    const char* varName,
    const string& value
)
{
    adios_define_attribute
    (
        groupID_,
        attrName,
        varName,
        adios_string,
        value.c_str(),
        NULL
    );
}


void Foam::adiosWrite::defineAttribute
(
    const char* attrName,
    const fileName& varName,
    const string& value
)
{
    defineAttribute(attrName, varName.c_str(), value);
}


void Foam::adiosWrite::defineIntAttribute
(
    const char* attrName,
    const char* varName,
    const int value
)
{
    int intval = value;

    adios_define_attribute_byvalue
    (
        groupID_,
        attrName,
        varName,
        adios_integer,
        1,
        &intval
    );
}


void Foam::adiosWrite::defineIntAttribute
(
    const char* attrName,
    const fileName& varName,
    const int value
)
{
    defineIntAttribute(attrName, varName.c_str(), value);
}


bool Foam::adiosWrite::defineListAttribute
(
    const char* attrName,
    const fileName& varName,
    const UList<int>& list
)
{
    if (list.empty())
    {
        return false;
    }

    adios_define_attribute_byvalue
    (
        groupID_,
        attrName,
        varName.c_str(),
        adios_integer,
        list.size(),
        const_cast<int*>(list.cdata())  // adios has "void*" instead of "const void*"
    );

    return true;
}


bool Foam::adiosWrite::defineListAttribute
(
    const char* attrName,
    const fileName& varName,
    const wordList& list
)
{
    CStringList cstrings;
    if (!list.empty())
    {
        cstrings.reset(list);
    }

    if (cstrings.argc())
    {
        adios_define_attribute_byvalue
        (
            groupID_,
            attrName,
            varName.c_str(),
            adios_string_array,
            cstrings.argc(),
            cstrings.argv()
        );

        return true;
    }

    return false;
}


bool Foam::adiosWrite::defineListAttribute
(
    const char* attrName,
    const fileName& varName,
    const stringList& list
)
{
    CStringList cstrings;
    if (!list.empty())
    {
        cstrings.reset(list);
    }

    if (cstrings.argc())
    {
        adios_define_attribute_byvalue
        (
            groupID_,
            attrName,
            varName.c_str(),
            adios_string_array,
            cstrings.argc(),
            cstrings.argv()
        );

        return true;
    }

    return false;
}


void Foam::adiosWrite::writeVariable
(
    const char* name,
    const void* value
)
{
    adios_write(fileID_, name, value);
}


void Foam::adiosWrite::writeVariable
(
    const fileName& name,
    const void* value
)
{
    writeVariable(name.c_str(), value);
}


void Foam::adiosWrite::writeIntVariable
(
    const char* name,
    const int value
)
{
    adios_write(fileID_, name, &value);
}


void Foam::adiosWrite::writeIntVariable
(
    const fileName& name,
    const int value
)
{
    writeIntVariable(name.c_str(), value);
}


void Foam::adiosWrite::writeScalarVariable
(
    const char* name,
    const double value
)
{
    adios_write(fileID_, name, &value);
}


void Foam::adiosWrite::writeScalarVariable
(
    const fileName& name,
    const double value
)
{
    writeScalarVariable(name.c_str(), value);
}


// ************************************************************************* //
