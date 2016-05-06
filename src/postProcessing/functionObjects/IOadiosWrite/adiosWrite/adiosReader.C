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

    for (int varI=0; varI < file->nvars; ++varI)
    {
        maxLen = Foam::max(maxLen, sizeOf(file->var_namelist[varI]));
    }

    if (verbose)
    {
        Info<< "max-length: " << maxLen << endl;
    }
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
    const fileName& datasetName,
    bool verbose
)
{
    return sizeOf(datasetName.c_str(), verbose);
}


bool Foam::adiosReader::helper::getDataSet
(
    const fileName& datasetName
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
    const fileName& datasetName,
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

    maxLen = 0;
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
