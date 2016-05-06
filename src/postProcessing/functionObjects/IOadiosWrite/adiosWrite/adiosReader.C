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


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::adiosReader::helper::scan(bool verbose)
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
        const char* attrName = file->attr_namelist[i];
        attributes.insert(attrName, i);

        if (verbose)
        {
            Info<< "attribute: " << attrName << endl;
        }
    }

    // variable:
    // hash as (name => nbytes on this process)
    for (int i=0; i < file->nvars; ++i)
    {
        const char* varName = file->var_namelist[i];

        size_t varsize = sizeOf(varName);
        variables.insert(varName, varsize);

        maxLen = Foam::max(maxLen, varsize);

        if (verbose)
        {
            Info<< "variable: " << varName << endl;
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
         && names.empty()
        )
        {
            cloudNames_.insert(regName, names);
        }
    }

    // for diagnostics: expect these type of entries
    // regionName/cloudName/nParticle
    // regionName/cloudName/size
    // regionName/cloudName/names
    // regionName/cloudName/types
    // regionName/cloudName/offset
    // regionName/cloudName/byte-size
    forAll(regionNames_, regI)
    {
        const word& regName = regionNames_[regI];
        if (cloudNames_.found(regName))
        {
            const wordList& cloudNames = cloudNames_[regName];

            forAll(cloudNames, cloudI)
            {
                fileName varPath = regName/"cloud"/cloudNames[cloudI];
                label nTotal = getIntAttribute(varPath/"nParticle");
                label nBytes = getIntAttribute(varPath/"size");

                wordList fragName = getStringListAttribute<word>(varPath/"names");
                wordList fragType = getStringListAttribute<word>(varPath/"types");

                labelList fragOff  = getIntListAttribute(varPath/"offset");
                labelList fragSize = getIntListAttribute(varPath/"byte-size");
            }

            // cloudNames_.insert(regName, names);
        }
    }

    if (verbose)
    {
        Info<< "max-length: " << maxLen << endl;
    }
}


Foam::HashTable<Foam::adiosReader::fieldInfo>
Foam::adiosReader::helper::getFieldInfo(const word& regName) const
{
    const string startsWith = regName / "field";

    typedef HashTable<size_t, fileName> VarContainer;

    HashTable<fieldInfo> table;

    forAllConstIter(VarContainer, variables, iter)
    {
        const fileName& varName = iter.key();
        size_t bytes = iter();

        if (varName.count('/') == 2 && varName.path() == startsWith)
        {
            // matches regName/field/xxx

            table.insert
            (
                varName.name(),
                fieldInfo
                (
                    varName,
                    bytes,
                    getStringAttribute(varName/"class")
                )
            );
        }
    }

    return table;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adiosReader::helper::helper(const DynamicCharList& buf)
:
    buffer(const_cast<DynamicCharList&>(buf)),
    file(NULL),
    selection(NULL),
    maxLen(0)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::adiosReader::helper::~helper()
{
    close();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


bool Foam::adiosReader::helper::open(const fileName& bpFile, MPI_Comm comm)
{
    // close anything already open
    close();

    Info<< " Read data of step " << bpFile << endl;

    file = adios_read_open_file
    (
        bpFile.c_str(),
        ADIOS_READ_METHOD_BP,
        comm
    );

    if (!file)
    {
        return false;
    }

    Info<< "found num-vars: " << file->nvars << endl;

    return file;
}


void Foam::adiosReader::helper::select(ADIOS_SELECTION *sel)
{
    if (selection)
    {
        adios_selection_delete(selection);
    }

    selection = sel;
}


size_t Foam::adiosReader::helper::sizeOf
(
    const char* datasetName,
    bool verbose
)
{
    size_t bytes = 0;

    if (file)
    {
        // TODO? restrict sizing to current processor!
        ADIOS_VARINFO *varInfo = adios_inq_var(file, datasetName);

        if (!varInfo)
        {
            WarningInFunction
                << "Error reading variable information " << datasetName
                << " from adios file: "
                << adios_errmsg() << endl;

            return 0;
        }

        if (varInfo->type == adios_string)
        {
            WarningInFunction
                << "Reading sizes for adios_string variables incomplete: " << datasetName
                << " from adios file: "
                << endl;

            return 0;
        }

        if (varInfo->ndim > 0)
        {
            int nblocks = varInfo->sum_nblocks;

            // Pout<< " variable=" << datasetName
            //     << " dims: " << varInfo->ndim << " nblocks:"  << nblocks << endl;

            // get block-decomposition
            int err = adios_inq_var_blockinfo(file, varInfo);
            if (err)
            {
                WarningInFunction
                    << "Error reading blockinfo for dataset " << datasetName
                    << " from adios file: "
                    << adios_errmsg() << endl;
            }
            else
            {
                ADIOS_VARBLOCK *bp = varInfo->blockinfo;

                for (int blockI=0; blockI < nblocks; ++blockI)
                {
                    // Pout<< datasetName
                    //     << " size:"<< varInfo->dims[0]
                    //     << " block[" << blockI
                    //     << "] start:" << *(bp->start)
                    //     << " count:" << *(bp->count)
                    //     << " process:" << bp->process_id
                    //     << " time-index:" << bp->time_index
                    //     << endl;

                    if (Foam::label(bp->process_id) == Pstream::myProcNo())
                    {
                        bytes += *(bp->count);
                    }

                    ++bp;
                }
            }

            // fallback: only consider 1-D storage:
            if (!bytes)
            {
                bytes = varInfo->dims[0];

                // for (int dimI=1; dimI < varInfo->ndim; ++dimI)
                // {
                //     ... varInfo->dims[dimI];
                // }
            }
        }
        else
        {
            // Pout<< " variable=" << datasetName << " is scalar" << endl;
            bytes = 1; // scalar value
        }

        bytes *= adios_type_size(varInfo->type, NULL);

        if (verbose)
        {
            Pout<< " variable=" << datasetName
                << " nbytes=" << bytes
                << " type=" << adios_type_to_string(varInfo->type)
                << endl;
        }

        // free ADIOS_VARINFO and any ADIOS_VARBLOCK(s)
        adios_free_varinfo(varInfo);
    }

    return bytes;
}


size_t Foam::adiosReader::helper::sizeOf
(
    const string& datasetName,
    bool verbose
)
{
    return sizeOf(datasetName.c_str(), verbose);
}


bool Foam::adiosReader::helper::getDataSet
(
    const string& datasetName
)
{
    // is.name() = Foam::name(Pstream::myProcNo()) / datasetName;
    size_t nbytes = sizeOf(datasetName, true);
    buffer.reserve(nbytes);

    bool ok = getDataSet(datasetName, buffer.data());
    if (ok)
    {
        buffer.setSize(nbytes);
    }
    else
    {
        buffer.setSize(0);
    }

    return ok;
}


bool Foam::adiosReader::helper::getDataSet
(
    const string& datasetName,
    void* data
)
{
    Pout<<"read data-set " << datasetName << endl;

    int err = adios_schedule_read
    (
        file,
        selection,
        datasetName.c_str(),
        0, 1,
        data
    );

    if (err)
    {
        WarningInFunction
            << "Error reading dataset " << datasetName
            << " from adios file: "
            << adios_errmsg() << endl;
    }
    else
    {
        err = adios_perform_reads(file, 1); // blocking
    }

    return !err;
}


void Foam::adiosReader::helper::close()
{
    if (selection)
    {
        adios_selection_delete(selection);
        selection = NULL;
    }

    if (file)
    {
        adios_read_close(file);
        file = NULL;
    }

    attributes.clear();
    variables.clear();

    regionNames_.clear();
    cloudNames_.clear();

    maxLen = 0;
}


bool Foam::adiosReader::helper::readIntAttributeIfPresent
(
    const string& attrName,
    label &value
) const
{
    enum ADIOS_DATATYPES type = adios_unknown;
    int  size = 0;
    void *data = 0;

    bool ok =
    (
        attributes.found(attrName)
     && 0 ==
        adios_get_attr_byid
        (
            file,
            attributes[attrName],
            &type,
            &size,
            &data
        )
    );

    if (ok)
    {
        if (type == adios_integer)
        {
            size /= sizeof(int);
            ok = (size == 1);

            if (ok)
            {
                value = *(reinterpret_cast<int*>(data));
            }
        }
        else if (type == adios_unsigned_integer)
        {
            size /= sizeof(unsigned int);
            ok = (size == 1);

            if (ok)
            {
                value = *(reinterpret_cast<unsigned int*>(data));
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
                    << "too many elements for attribute: " << attrName << nl
                    << "  expecting 1 int/unsigned, found " << size
                    << exit(FatalIOError);
            }
            else
            {
                FatalErrorInFunction
                    << "incorrect type for attribute: " << attrName << nl
                    << "  expecting int/unsigned, found "
                    << adios_type_to_string(type)
                    << exit(FatalIOError);
            }
        }
    }

    return ok;
}


bool Foam::adiosReader::helper::readIntListAttributeIfPresent
(
    const string& attrName,
    List<label>& lst
) const
{
    enum ADIOS_DATATYPES type = adios_unknown;
    int  size = 0;
    void *data = 0;

    bool ok =
    (
        attributes.found(attrName)
     && 0 ==
        adios_get_attr_byid
        (
            file,
            attributes[attrName],
            &type,
            &size,
            &data
        )
    );

    if (ok)
    {
        if (type == adios_integer)
        {
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
        else if (type == adios_unsigned_integer)
        {
            size /= sizeof(unsigned int);
            ok = (size > 1);

            if (ok)
            {
                unsigned int* ptr = reinterpret_cast<unsigned int*>(data);

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
                    << "too few elements for attribute: " << attrName << nl
                    << "  expecting multiple int/unsigned values, found 1"
                    << exit(FatalIOError);
            }
            else
            {
                FatalErrorInFunction
                    << "incorrect type for attribute: " << attrName << nl
                    << "  expecting int/unsigned, found "
                    << adios_type_to_string(type)
                    << exit(FatalIOError);
            }
        }
    }

    return ok;
}


bool Foam::adiosReader::helper::readStringAttributeIfPresent
(
    const string& attrName,
    string &value
) const
{
    enum ADIOS_DATATYPES type = adios_unknown;
    int  size = 0;
    void *data = 0;

    bool ok =
    (
        attributes.found(attrName)
     && 0 ==
        adios_get_attr_byid
        (
            file,
            attributes[attrName],
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
                << "incorrect type for attribute: " << attrName << nl
                << "  expecting string, found "
                << adios_type_to_string(type)
                << exit(FatalIOError);
        }
    }

    return ok;
}


Foam::label Foam::adiosReader::helper::getIntAttribute
(
    const string& attrName
) const
{
    label value;

    if (!readIntAttributeIfPresent(attrName, value))
    {
        FatalErrorInFunction
            << "integer attribute missing: " << attrName
            << exit(FatalIOError);
    }

    return value;
}


Foam::List<Foam::label> Foam::adiosReader::helper::getIntListAttribute
(
    const string& attrName
) const
{
    List<label> value;

    if (!readIntListAttributeIfPresent(attrName, value))
    {
        FatalErrorInFunction
            << "int-list attribute missing: " << attrName
            << exit(FatalIOError);
    }

    return value;
}


Foam::string Foam::adiosReader::helper::getStringAttribute
(
    const string& attrName
) const
{
    string value;

    if (!readStringAttributeIfPresent(attrName, value))
    {
        FatalErrorInFunction
            << "string attribute missing: " << attrName
            << exit(FatalIOError);
    }

    return value;
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
