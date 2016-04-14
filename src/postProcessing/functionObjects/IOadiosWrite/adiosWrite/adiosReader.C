/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd
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
#include "dictionary.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::adiosReader::helper::scan(bool verbose)
{
    maxLen = 0;

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
        Info<< "adios-file has " << file->nvars << " variables" << endl;
    }


    // TODO? restrict sizing to current processor!
    for (int varI=0; varI < file->nvars; ++varI)
    {
        const char * varName = file->var_namelist[varI];
        ADIOS_VARINFO *varInfo = adios_inq_var(file, varName);

        if (!varInfo)
        {
            WarningInFunction
                << "Error reading variable information " << varName
                << " from adios file: "
                << adios_errmsg() << endl;
            continue;
        }

        size_t bytes = 1; // fallback to scalar
        if (varInfo->ndim > 0)
        {
            // only consider 1-D storage:
            bytes = varInfo->dims[0];
            // for (int dimI=1; dimI < varInfo->ndim; ++dimI)
            // {
            //     ... varInfo->dims[dimI];
            // }
        }

        bytes *= adios_type_size(varInfo->type, const_cast<char *>(""));

        if (verbose)
        {
            Info<< " variable=" << varName
                << " nbytes=" << bytes
                << " type=" << adios_type_to_string(varInfo->type)
                << endl;
        }

        maxLen = Foam::max(maxLen, bytes);

        // free ADIOS_VARINFO
        adios_free_varinfo(varInfo);
    }

    if (verbose)
    {
        Info<< "max-length: " << maxLen << endl;
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adiosReader::helper::helper()
:
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


bool Foam::adiosReader::helper::getDataSet
(
    const fileName& datasetName,
    void* data
)
{
    Info<<"read data-set " << datasetName << endl;

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


bool Foam::adiosReader::helper::getDataSet
(
    const fileName& datasetName,
    IStringStreamBuf& is
)
{
    is.rewind();
    return getDataSet(datasetName, is.data());
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

    maxLen = 0;
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
