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

#include "adiosReader.H"
#include "adiosCore.H"

#include "dictionary.H"
#include "IOstreams.H"
#include "Pstream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::adiosReader::scan(bool verbose)
{
    maxLen = 0;
    variables.clear();
    attributes.clear();

    regionNames_.clear();
    cloudNames_.clear();

    if (!file)
    {
        if (verbose)
        {
            Info<< "adios-file not open" << endl;
        }
        return;
    }

    if (verbose)
    {
        Info<< "adios-file has " << file->nvars << " variables and "
            << file->nattrs << " attributes" << endl;
    }


    // attributes:
    // hash as (name => lookup index)
    for (int i=0; i < file->nattrs; ++i)
    {
        const char* name = file->attr_namelist[i];
        attributes.insert(name, i);

        if (verbose)
        {
            Info<< "attribute: " << name << endl;
        }
    }

    // variable:
    // hash as (name => nbytes on this process)
    for (int i=0; i < file->nvars; ++i)
    {
        VarInfo vinfo(file, file->var_namelist[i]);

        variables.insert(vinfo.name(), vinfo);

        maxLen = Foam::max(maxLen, vinfo.sizeOf());

        if (verbose)
        {
            Info<< "variable: " << vinfo.name() << endl;
        }
    }


    // mandatory: /openfoam/nRegions - but can also just rely on region names

    // mandatory: /openfoam/regions
    regionNames_ = getStringListAttribute<word>
    (
        adiosCore::foamAttribute/"regions"
    );

    // optional: regionName/nClouds, regionName/clouds
    forAll(regionNames_, regI)
    {
        const word& regName = regionNames_[regI];

        wordList names;
        if
        (
            readStringListAttributeIfPresent(regName/"clouds", names)
         && !names.empty()
        )
        {
            cloudNames_.insert(regName, names);
        }
    }

    if (verbose)
    {
        Info<< "max-length: " << maxLen << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adiosReader::adiosReader()
:
    file(0),
    selection(0),
    maxLen(0),
    regionNames_(),
    cloudNames_()
{}


Foam::adiosReader::adiosReader(const fileName& bpFile, MPI_Comm comm)
:
    file(0),
    selection(0),
    maxLen(0),
    regionNames_(),
    cloudNames_()
{
    open(bpFile, comm);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::adiosReader::~adiosReader()
{
    close();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


bool Foam::adiosReader::open(const fileName& bpFile, MPI_Comm comm)
{
    // close anything already open
    close();

    // Info<< " Read data of step " << bpFile << endl;

    file = adios_read_open_file
    (
        bpFile.c_str(),
        ADIOS_READ_METHOD_BP,
        comm
    );

    if (file)
    {
        select(adios_selection_writeblock(Pstream::myProcNo()));
        scan(false);

        // Info<< "found num-vars: " << file->nvars << endl;
    }
    else
    {
        return false;
    }

    return isGood();
}


void Foam::adiosReader::select(ADIOS_SELECTION *sel)
{
    if (selection)
    {
        adios_selection_delete(selection);
    }

    selection = sel;
}


void Foam::adiosReader::close()
{
    if (selection)
    {
        adios_selection_delete(selection);
        selection = 0;
    }

    if (file)
    {
        adios_read_close(file);
        file = 0;
    }

    attributes.clear();
    variables.clear();

    regionNames_.clear();
    cloudNames_.clear();

    maxLen = 0;
}


bool Foam::adiosReader::isGood() const
{
    return file;
}


Foam::HashTable<Foam::adiosReader::fieldInfo>
Foam::adiosReader::getFieldInfo(const word& regName) const
{
    const string startsWith = adiosCore::fieldPath(regName);

    HashTable<fieldInfo> table;

    forAllConstIter(VarContainer, variables, iter)
    {
        const fileName& varName = iter.key();
        const VarInfo&  varInfo = iter();

        if (varName.count('/') == 2 && varName.path() == startsWith)
        {
            // matches regName/field/xxx

            fieldInfo info
            (
                varName,
                varInfo.sizeOf(),
                getStringAttribute(varName/"class")
            );

            table.insert(varName.name(), info);
        }
    }

    return table;
}


#ifdef HAS_CLOUD_SUPPORT
Foam::HashTable<Foam::adiosReader::cloudInfo>
Foam::adiosReader::getCloudInfo(const word& regName) const
{
    const string startsWith = adiosCore::cloudPath(regName);

    HashTable<cloudInfo> table;

    forAllConstIter(VarContainer, variables, iter)
    {
        const fileName& varName = iter.key();
        const VarInfo&  varInfo = iter();

        if (varName.count('/') == 2 && varName.path() == startsWith)
        {
            // matches regName/cloud/xxx
            table.insert(varName.name(), cloudInfo(varInfo, *this));
        }
    }

    return table;
}
#endif

bool Foam::adiosReader::readIntAttributeIfPresent
(
    const string& name,
    label& value
) const
{
    enum ADIOS_DATATYPES type = adios_unknown;
    int  size = 0;
    void *data = 0;

    bool ok =
    (
        hasAttribute(name)
     && 0 ==
        adios_get_attr_byid
        (
            file,
            attributes[name],
            &type,
            &size,
            &data
        )
    );

    if (ok)
    {
        if (type == adios_integer || type == adios_unsigned_integer)
        {
            // we don't distinguish between signed/unsigned
            // mostly just use signed anyhow

            size /= sizeof(int);
            ok = (size == 1);

            if (ok)
            {
                value = *(reinterpret_cast<int*>(data));
            }
        }
        else
        {
            size = 0;
            ok = false;
        }

        if (data)
        {
            free(data);
        }

        if (!ok)
        {
            if (size)
            {
                FatalErrorInFunction
                    << "too many elements for attribute: " << name << nl
                    << "  expecting 1 int/unsigned, found " << size
                    << exit(FatalIOError);
            }
            else
            {
                FatalErrorInFunction
                    << "incorrect type for attribute: " << name << nl
                    << "  expecting int/unsigned, found "
                    << adios_type_to_string(type)
                    << exit(FatalIOError);
            }
        }
    }

    return ok;
}


bool Foam::adiosReader::readScalarAttributeIfPresent
(
    const string& name,
    double& value
) const
{
    enum ADIOS_DATATYPES type = adios_unknown;
    int  size = 0;
    void *data = 0;

    bool ok =
    (
        hasAttribute(name)
     && 0 ==
        adios_get_attr_byid
        (
            file,
            attributes[name],
            &type,
            &size,
            &data
        )
    );

    if (ok)
    {
        if (type == adios_double)
        {
            size /= sizeof(double);
            ok = (size == 1);

            if (ok)
            {
                value = *(reinterpret_cast<int*>(data));
            }
        }
        else
        {
            size = 0;
            ok = false;
        }

        if (data)
        {
            free(data);
        }

        if (!ok)
        {
            if (size)
            {
                FatalErrorInFunction
                    << "too many elements for attribute: " << name << nl
                    << "  expecting 1 double, found " << size
                    << exit(FatalIOError);
            }
            else
            {
                FatalErrorInFunction
                    << "incorrect type for attribute: " << name << nl
                    << "  expecting double, found "
                    << adios_type_to_string(type)
                    << exit(FatalIOError);
            }
        }
    }

    return ok;
}


bool Foam::adiosReader::readScalarListAttributeIfPresent
(
    const string& name,
    List<double>& lst
) const
{
    enum ADIOS_DATATYPES type = adios_unknown;
    int  size = 0;
    void *data = 0;

    bool ok =
    (
        hasAttribute(name)
     && 0 ==
        adios_get_attr_byid
        (
            file,
            attributes[name],
            &type,
            &size,
            &data
        )
    );

    if (ok)
    {
        if (type == adios_double)
        {
            size /= sizeof(double);
            ok = (size > 1);

            if (ok)
            {
                double *ptr = reinterpret_cast<double*>(data);

                lst.setSize(size);
                for (int i=0; i < size; ++i)
                {
                    lst[i] = ptr[i];
                }
            }
        }
        else
        {
            size = 0;
            ok = false;
        }

        if (data)
        {
            free(data);
        }

        if (!ok)
        {
            if (size == 1)
            {
                FatalErrorInFunction
                    << "too few elements for attribute: " << name << nl
                    << "  expecting multiple double values, found 1"
                    << exit(FatalIOError);
            }
            else
            {
                FatalErrorInFunction
                    << "incorrect type for attribute: " << name << nl
                    << "  expecting double, found "
                    << adios_type_to_string(type)
                    << exit(FatalIOError);
            }
        }
    }

    return ok;
}


bool Foam::adiosReader::readIntListAttributeIfPresent
(
    const string& name,
    List<label>& lst
) const
{
    enum ADIOS_DATATYPES type = adios_unknown;
    int  size = 0;
    void *data = 0;

    bool ok =
    (
        hasAttribute(name)
     && 0 ==
        adios_get_attr_byid
        (
            file,
            attributes[name],
            &type,
            &size,
            &data
        )
    );

    if (ok)
    {
        if (type == adios_integer || adios_unsigned_integer)
        {
            // we don't distinguish between signed/unsigned
            // mostly just use signed anyhow

            size /= sizeof(int);
            ok = (size > 1);

            if (ok)
            {
                int *ptr = reinterpret_cast<int*>(data);

                lst.setSize(size);
                for (int i=0; i < size; ++i)
                {
                    lst[i] = ptr[i];
                }
            }
        }
        else
        {
            size = 0;
            ok = false;
        }

        if (data)
        {
            free(data);
        }

        if (!ok)
        {
            if (size == 1)
            {
                FatalErrorInFunction
                    << "too few elements for attribute: " << name << nl
                    << "  expecting multiple int/unsigned values, found 1"
                    << exit(FatalIOError);
            }
            else
            {
                FatalErrorInFunction
                    << "incorrect type for attribute: " << name << nl
                    << "  expecting int/unsigned, found "
                    << adios_type_to_string(type)
                    << exit(FatalIOError);
            }
        }
    }

    return ok;
}


bool Foam::adiosReader::readStringAttributeIfPresent
(
    const string& name,
    string& value
) const
{
    enum ADIOS_DATATYPES type = adios_unknown;
    int  size = 0;
    void *data = 0;

    bool ok =
    (
        hasAttribute(name)
     && 0 ==
        adios_get_attr_byid
        (
            file,
            attributes[name],
            &type,
            &size,
            &data
        )
    );

    if (ok)
    {
        if (type == adios_string)
        {
            value = reinterpret_cast<char*>(data);
        }
        else
        {
            ok = false;
        }

        if (data)
        {
            free(data);
        }

        if (!ok)
        {
            FatalErrorInFunction
                << "incorrect type for attribute: " << name << nl
                << "  expecting string, found "
                << adios_type_to_string(type)
                << exit(FatalIOError);
        }
    }

    return ok;
}


bool Foam::adiosReader::readIntVariableIfPresent
(
    const string& name,
    label& value
) const
{
    bool ok = hasVariable(name);
    int err = 0;
    int data = 0;

    if (ok)
    {
        const VarInfo& vinfo = variables[name];
        int size = vinfo.nElem();

        if
        (
            vinfo.dataType() == adios_integer
         || vinfo.dataType() == adios_unsigned_integer
        )
        {
            // we don't distinguish between signed/unsigned
            // mostly just use signed anyhow

            ok = (size == 1);

            if (ok)
            {
                err = adios_schedule_read
                (
                    file,
                    selection,
                    name.c_str(),
                    0, 1,
                    &data
                );

                if (!err)
                {
                    err = adios_perform_reads(file, 1); // blocking
                }

                if (!err)
                {
                    value = data;
                }
                else
                {
                    FatalErrorInFunction
                        << "error reading adios variable: " << name << nl
                        << exit(FatalIOError);
                }
            }
        }
        else
        {
            size = 0;
            ok = false;
        }

        if (!ok)
        {
            if (size)
            {
                FatalErrorInFunction
                    << "too many elements for variable: " << name << nl
                    << "  expecting 1 int/unsigned, found " << size
                    << exit(FatalIOError);
            }
            else
            {
                FatalErrorInFunction
                    << "incorrect type for variable: " << name << nl
                    << "  expecting int/unsigned, found "
                    << adios_type_to_string(vinfo.dataType())
                    << exit(FatalIOError);
            }
        }
    }

    return ok;
}


bool Foam::adiosReader::readScalarVariableIfPresent
(
    const string& name,
    scalar& value
) const
{
    bool ok = hasVariable(name);
    int err = 0;

    if (ok)
    {
        const VarInfo& vinfo = variables[name];

        int size = vinfo.nElem();
        if (vinfo.dataType() == adios_real)
        {
            ok = (size == 1);

            if (ok)
            {
                float data = 0;

                err = adios_schedule_read
                (
                    file,
                    selection,
                    name.c_str(),
                    0, 1,
                    &data
                );

                if (!err)
                {
                    err = adios_perform_reads(file, 1); // blocking
                }

                if (!err)
                {
                    value = data;
                }
                else
                {
                    FatalErrorInFunction
                        << "error reading adios variable: " << name << nl
                        << exit(FatalIOError);
                }
            }
        }
        else if (vinfo.dataType() == adios_double)
        {
            ok = (size == 1);

            if (ok)
            {
                double data = 0;

                err = adios_schedule_read
                (
                    file,
                    selection,
                    name.c_str(),
                    0, 1,
                    &data
                );

                if (!err)
                {
                    err = adios_perform_reads(file, 1); // blocking
                }

                if (!err)
                {
                    value = data;
                }
                else
                {
                    FatalErrorInFunction
                        << "error reading adios variable: " << name << nl
                        << exit(FatalIOError);
                }
            }
        }
        else
        {
            size = 0;
            ok = false;
        }

        if (!ok)
        {
            if (size)
            {
                FatalErrorInFunction
                    << "too many elements for variable: " << name << nl
                    << "  expecting 1 float/double, found " << size
                    << exit(FatalIOError);
            }
            else
            {
                FatalErrorInFunction
                    << "incorrect type for variable: " << name << nl
                    << "  expecting  float/double, found "
                    << adios_type_to_string(vinfo.dataType())
                    << exit(FatalIOError);
            }
        }
    }

    return ok;
}


bool Foam::adiosReader::getVariable
(
    const string& name,
    void* data
) const
{
    // Pout<<"read data-set " << name << endl;
    bool ok = hasVariable(name);
    if (ok)
    {
        int err = adios_schedule_read
        (
            file,
            selection,
            name.c_str(),
            0, 1,
            data
        );

        if (!err)
        {
            err = adios_perform_reads(file, 1); // blocking
        }

        if (err)
        {
            FatalErrorInFunction
                << "Error reading adios variable " << name
                << adios_errmsg() << endl;
        }

        ok = !err;
    }
    else
    {
        FatalErrorInFunction
            << "missing adios variable: " << name
            << exit(FatalIOError);
    }

    return ok;
}


// ************************************************************************* //
