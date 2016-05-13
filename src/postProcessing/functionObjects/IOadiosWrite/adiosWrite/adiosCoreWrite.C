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

#include "adiosCoreWrite.H"
#include "adiosTime.H"
#include "foamVersion.H"
#include "OSspecific.H"

// some internal pre-processor stringifications
#undef STRINGIFY
#undef TO_STRING
#define STRINGIFY(x) #x
#define TO_STRING(x) STRINGIFY(x)


// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

void Foam::adiosCoreWrite::read(const dictionary& dict)
{
    // Default values prior to lookup in dictionary
    readMethod_  = "BP";
    writeMethod_ = "MPI";
    writeParams_ = "";

    dict.readIfPresent("readMethod",  readMethod_);
    dict.readIfPresent("adiosMethod", writeMethod_);  // old name
    dict.readIfPresent("writeMethod", writeMethod_);
    dict.readIfPresent("methodparams", writeParams_); // old name
    dict.readIfPresent("writeOptions", writeParams_);

    Info<< "  ADIOS writeMethod: " << writeMethod_ << nl
        << "        writeParams: " << writeParams_ << nl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adiosCoreWrite::adiosCoreWrite
(
    const word& groupName,
    const dictionary& dict
)
:
    adiosCore(groupName),
    groupID_(0),
    fileID_(0),
    readMethod_(),
    writeMethod_(),
    writeParams_(),
    comm_(MPI_COMM_NULL),
    vars_()
{
    if (!Pstream::parRun())
    {
        // serial - must do MPI_Init() ourselves
        MPI_Init(NULL, NULL);  // NULL args are OK for openmpi, mpich-2
    }

    MPI_Comm_dup(MPI_COMM_WORLD, &comm_);

    // Initialize ADIOS
    adios_init_noxml(comm_);

    // Create a group to hold all variable definitions
    adios_declare_group(&groupID_, name().c_str(), "", adios_flag_yes);

    // Read dictionary
    read(dict);

    // Set the actual output method here for the ADIOS group
    adios_select_method
    (
        groupID_,
        writeMethod_.c_str(),
        writeParams_.c_str(),
        ""  // base-path (unused) needs empty string, not a NULL pointer
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::adiosCoreWrite::~adiosCoreWrite()
{
    adios_free_group(groupID_); // not necessary but nice to cleanup
    adios_finalize(Pstream::myProcNo());

    MPI_Comm_free(&comm_);

    if (!Pstream::parRun())
    {
        // serial - must do MPI_Finalize() ourselves
        MPI_Finalize();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::adiosCoreWrite::open(const fileName& dataFile)
{
    fileName path = dataFile.path();
    if (!isDir(path))
    {
        mkDir(path);
    }

    char mode[] = "w"; // create filename
    int err = adios_open
    (
        &fileID_,
        name().c_str(),     // group-name
        dataFile.c_str(),
        mode,
        comm_
    );

    // Info<< "adiosCoreWrite::open (" << dataFile << ")" << endl;
    if (err)
    {
        // OR Fatal?
        WarningInFunction
            << "File " << dataFile
            << " could not be created: " << nl
            << adios_get_last_errmsg();

        fileID_ = 0;
    }

    return fileID_ != 0;
}


void Foam::adiosCoreWrite::close()
{
    MPI_Barrier(MPI_COMM_WORLD);
    // Info<< "adiosCoreWrite::close" << nl << endl;

    // Close the file
    if (fileID_)
    {
        adios_close(fileID_);
        fileID_ = 0;
    }
}


void Foam::adiosCoreWrite::reset()
{
    // In ADIOS we need to remove all variable definitions in order
    // to make a new list of definitions in case the mesh changes
    adios_delete_vardefs(groupID_);

    // also cleanup all old attributes
    adios_delete_attrdefs(groupID_);

    vars_.clear(); // cleanup name => id mappings
}


size_t Foam::adiosCoreWrite::maxSize() const
{
    size_t maxLen = 0;

    forAllConstIter(VarHashType, vars_, iter)
    {
        maxLen = Foam::max(maxLen, adios_expected_var_size(iter()));
    }

    return maxLen;
}


size_t Foam::adiosCoreWrite::outputSize() const
{
    size_t total = 0;

    forAllConstIter(VarHashType, vars_, iter)
    {
        total += adios_expected_var_size(iter());
    }

    return total;
}


void Foam::adiosCoreWrite::defineBaseAttributes()
{
    // OpenFOAM build information (also contains version)
    defineAttribute
    (
        "version",
        adiosCore::foamAttribute,
        Foam::FOAMbuild
    );

    // OpenFOAM platform tag (WM_ARCH + WM_COMPILER)
    defineAttribute
    (
        "platform",
        adiosCore::foamAttribute,
#ifdef FOAM_PLATFORM
        TO_STRING(FOAM_PLATFORM)
#else
        ""
# error "FOAM_PLATFORM not defined"
#endif
    );

    // OpenFOAM label size (32|64)
    defineIntAttribute
    (
        "label-size",
        adiosCore::foamAttribute,
        adiosTraits<label>::nBits
    );

    // OpenFOAM scalar type (single-precision?)
    defineAttribute
    (
        "precision",
        adiosCore::foamAttribute,
        adiosTraits<scalar>::precisionName
    );

    // other general information
    defineIntAttribute("nProcs", adiosCore::foamAttribute, Pstream::nProcs());
}


void Foam::adiosCoreWrite::defineTimeAttributes(const TimeState& t)
{
    defineIntAttribute
    (
        adiosTime::attr[adiosTime::INDEX],
        adiosCore::timeAttribute,
        t.timeIndex()
    );
    defineScalarAttribute
    (
        adiosTime::attr[adiosTime::VALUE],
        adiosCore::timeAttribute,
        t.timeOutputValue()
    );
    defineScalarAttribute
    (
        adiosTime::attr[adiosTime::DT],
        adiosCore::timeAttribute,
        t.deltaTValue()
    );
    defineScalarAttribute
    (
        adiosTime::attr[adiosTime::DT0],
        adiosCore::timeAttribute,
        t.deltaT0Value()
    );
}


void Foam::adiosCoreWrite::defineTime()
{
    defineIntVariable(adiosCore::timeAttribute/adiosTime::attr[adiosTime::INDEX]);
    defineScalarVariable(adiosCore::timeAttribute/adiosTime::attr[adiosTime::VALUE]);
    defineScalarVariable(adiosCore::timeAttribute/adiosTime::attr[adiosTime::DT]);
    defineScalarVariable(adiosCore::timeAttribute/adiosTime::attr[adiosTime::DT0]);
}


void Foam::adiosCoreWrite::writeTime(const TimeState& t)
{
    writeIntVariable
    (
        adiosCore::timeAttribute/adiosTime::attr[adiosTime::INDEX],
        t.timeIndex()
    );
    writeScalarVariable
    (
        adiosCore::timeAttribute/adiosTime::attr[adiosTime::VALUE],
        t.timeOutputValue()
    );
    writeScalarVariable
    (
        adiosCore::timeAttribute/adiosTime::attr[adiosTime::DT],
        t.deltaTValue()
    );
    writeScalarVariable
    (
        adiosCore::timeAttribute/adiosTime::attr[adiosTime::DT0],
        t.deltaT0Value()
    );
}


int64_t Foam::adiosCoreWrite::defineVariable
(
    const string& name,
    enum ADIOS_DATATYPES type,
    const string& local,
    const string& global,
    const string& offsets
)
{
    int64_t varid = adios_define_var
    (
        groupID_,
        name.c_str(),           // name
        NULL,                   // path (deprecated)
        type,                   // data-type
        local.c_str(),          // local dimensions
        global.c_str(),         // global dimensions
        offsets.c_str()         // local offsets
    );
    vars_.insert(name, varid);
    return varid;
}


int64_t Foam::adiosCoreWrite::defineVariable
(
    const char* name,
    enum ADIOS_DATATYPES type,
    size_t count
)
{
    int64_t varid = adios_define_var
    (
        groupID_,
        name,                           // name
        NULL,                           // path (deprecated)
        type,                           // data-type
        Foam::name(count).c_str(),      // local dimensions
        NULL,                           // global dimensions
        NULL                            // local offsets
    );
    vars_.insert(name, varid);
    return varid;
}


int64_t Foam::adiosCoreWrite::defineVariable
(
    const string& name,
    enum ADIOS_DATATYPES type,
    size_t count
)
{
    return defineVariable(name.c_str(), type, count);
}


int64_t Foam::adiosCoreWrite::defineIntVariable
(
    const char* name
)
{
    int64_t varid = adios_define_var
    (
        groupID_,
        name,                           // name
        NULL,                           // path (deprecated)
        adios_integer,                  // data-type
        "1",                            // local dimensions
        Foam::name(Pstream::nProcs()).c_str(),  // global 1D array of this info
        Foam::name(Pstream::myProcNo()).c_str() // offsets of this process into array
    );
    vars_.insert(name, varid);
    return varid;
}


int64_t Foam::adiosCoreWrite::defineIntVariable
(
    const string& name
)
{
    return defineIntVariable(name.c_str());
}


int64_t Foam::adiosCoreWrite::defineIntVariable
(
    const char* name,
    size_t count
)
{
    int64_t varid = adios_define_var
    (
        groupID_,
        name,                           // name
        NULL,                           // path (deprecated)
        adiosTraits<label>::adiosType,  // data-type
        Foam::name(count).c_str(),      // local dimensions
        NULL,                           // global dimensions
        NULL                            // local offsets
    );
    vars_.insert(name, varid);
    return varid;
}


int64_t Foam::adiosCoreWrite::defineIntVariable
(
    const string& name,
    size_t count
)
{
    return defineIntVariable(name.c_str(), count);
}


int64_t Foam::adiosCoreWrite::defineScalarVariable
(
    const char* name
)
{
    int64_t varid = adios_define_var
    (
        groupID_,
        name,                           // name
        NULL,                           // path (deprecated)
        adiosTraits<scalar>::adiosType, // data-type
        "1",                            // local dimensions
        Foam::name(Pstream::nProcs()).c_str(),  // global 1D array of this info
        Foam::name(Pstream::myProcNo()).c_str() // offsets of this process into array
    );
    vars_.insert(name, varid);
    return varid;
}


int64_t Foam::adiosCoreWrite::defineScalarVariable
(
    const string& name
)
{
    return defineScalarVariable(name.c_str());
}


int64_t Foam::adiosCoreWrite::defineScalarVariable
(
    const char* name,
    size_t count
)
{
    int64_t varid = adios_define_var
    (
        groupID_,
        name,                           // name
        NULL,                           // path (deprecated)
        adiosTraits<scalar>::adiosType, // data-type
        Foam::name(count).c_str(),      // local dimensions
        NULL,                           // global dimensions
        NULL                            // local offsets
    );
    vars_.insert(name, varid);
    return varid;
}


int64_t Foam::adiosCoreWrite::defineScalarVariable
(
    const string& name,
    size_t count
)
{
    return defineScalarVariable(name.c_str(), count);
}


int64_t Foam::adiosCoreWrite::defineStreamVariable
(
    const char* name,
    size_t count
)
{
    // use unsigned byte:
    // keeps people from thinking that the min/max statistics have any meaning

    int64_t varid = adios_define_var
    (
        groupID_,
        name,                           // name
        NULL,                           // path (deprecated)
        adios_unsigned_byte,            // data-type
        Foam::name(count).c_str(),      // local dimensions
        NULL,                           // global dimensions
        NULL                            // local offsets
    );
    vars_.insert(name, varid);
    return varid;
}


int64_t Foam::adiosCoreWrite::defineStreamVariable
(
    const string& name,
    size_t count
)
{
    return defineStreamVariable(name.c_str(), count);
}


int64_t Foam::adiosCoreWrite::defineVectorVariable
(
    const char* name
)
{
    const string global = Foam::name(Pstream::nProcs()) + ",3";
    const string offset = Foam::name(Pstream::myProcNo());

    int64_t varid = adios_define_var
    (
        groupID_,
        name,                           // name
        NULL,                           // path (deprecated)
        adiosTraits<scalar>::adiosType, // data-type
        "1,3",                          // local dimensions
        global.c_str(),                 // global dimensions
        offset.c_str()                  // offsets of this process into array
    );
    vars_.insert(name, varid);
    return varid;
}


int64_t Foam::adiosCoreWrite::defineVectorVariable
(
    const string& name
)
{
    return defineVectorVariable(name.c_str());
}


int64_t Foam::adiosCoreWrite::defineVectorVariable
(
    const char* name,
    size_t count
)
{
    const string local = Foam::name(count) + ",3";

    int64_t varid = adios_define_var
    (
        groupID_,
        name,                           // name
        NULL,                           // path (deprecated)
        adiosTraits<scalar>::adiosType, // data-type
        local.c_str(),                  // local dimensions
        NULL,                           // global dimensions
        NULL                            // local offsets
    );
    vars_.insert(name, varid);
    return varid;
}


int64_t Foam::adiosCoreWrite::defineVectorVariable
(
    const string& name,
    size_t count
)
{
    return defineVectorVariable(name.c_str(), count);
}


void Foam::adiosCoreWrite::defineAttribute
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


void Foam::adiosCoreWrite::defineAttribute
(
    const char* attrName,
    const string& varName,
    const string& value
)
{
    defineAttribute(attrName, varName.c_str(), value);
}


void Foam::adiosCoreWrite::defineIntAttribute
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


void Foam::adiosCoreWrite::defineIntAttribute
(
    const char* attrName,
    const string& varName,
    const int value
)
{
    defineIntAttribute(attrName, varName.c_str(), value);
}


void Foam::adiosCoreWrite::defineScalarAttribute
(
    const char* attrName,
    const char* varName,
    const double value
)
{
    double fltval = value;

    adios_define_attribute_byvalue
    (
        groupID_,
        attrName,
        varName,
        adios_double,
        1,
        &fltval
    );
}


void Foam::adiosCoreWrite::defineScalarAttribute
(
    const char* attrName,
    const string& varName,
    const double value
)
{
    defineScalarAttribute(attrName, varName.c_str(), value);
}


bool Foam::adiosCoreWrite::defineListAttribute
(
    const char* attrName,
    const string& varName,
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
        list.cdata()
    );

    return true;
}


bool Foam::adiosCoreWrite::defineListAttribute
(
    const char* attrName,
    const string& varName,
    const UList<double>& list
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
        adios_double,
        list.size(),
        list.cdata()
    );

    return true;
}


bool Foam::adiosCoreWrite::defineListAttribute
(
    const char* attrName,
    const string& varName,
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


bool Foam::adiosCoreWrite::defineListAttribute
(
    const char* attrName,
    const string& varName,
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


void Foam::adiosCoreWrite::writeVariable
(
    const char* name,
    const void* value
)
{
    // could also check that variable has been defined
    // vars_.found(name);

    if (fileID_)
    {
        adios_write(fileID_, name, value);
    }
    else
    {
        WarningInFunction
            << "Attempting to write adios variable \""
            << name << "\" without an open adios file"
            << endl;
    }
}


void Foam::adiosCoreWrite::writeVariable
(
    const string& name,
    const void* value
)
{
    writeVariable(name.c_str(), value);
}


void Foam::adiosCoreWrite::writeIntVariable
(
    const char* name,
    const int value
)
{
    // could also check that variable has been defined
    // vars_.found(name);
    if (fileID_)
    {
        adios_write(fileID_, name, &value);
    }
    else
    {
        WarningInFunction
            << "Attempting to write adios variable \""
            << name << "\" without an open adios file"
            << endl;
    }
}


void Foam::adiosCoreWrite::writeIntVariable
(
    const string& name,
    const int value
)
{
    writeIntVariable(name.c_str(), value);
}


void Foam::adiosCoreWrite::writeScalarVariable
(
    const char* name,
    const double value
)
{
    // could also check that variable has been defined
    // vars_.found(name);
    if (fileID_)
    {
        adios_write(fileID_, name, &value);
    }
    else
    {
        WarningInFunction
            << "Attempting to write adios variable \""
            << name << "\" without an open adios file"
            << endl;
    }
}


void Foam::adiosCoreWrite::writeScalarVariable
(
    const string& name,
    const double value
)
{
    writeScalarVariable(name.c_str(), value);
}


// ************************************************************************* //
