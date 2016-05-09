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

#include "IOstreams.H"
#include "Pstream.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::adiosReader::VarInfo::readInfo
(
    ADIOS_FILE *file,
    bool verbose
)
{
    nElem_  = 0;  // the number of elements (flattened dimensions)
    nBytes_ = 0;

    if (file)
    {
        // TODO? restrict sizing to current processor!
        ADIOS_VARINFO *vinfo = adios_inq_var(file, name_.c_str());

        if (!vinfo)
        {
            WarningInFunction
                << "Error reading variable information " << name_
                << " from adios file: "
                << adios_errmsg() << endl;

            return false;
        }

        if (vinfo->type == adios_string)
        {
            WarningInFunction
                << "Reading sizes for adios_string variables incomplete: "
                << name_ << endl;


            // free ADIOS_VARINFO and any ADIOS_VARBLOCK(s)
            adios_free_varinfo(vinfo);

            return false;
        }


        // could also assert (vinfo->nstep == 1)
        // since that is what we expect and support

        varid_ = vinfo->varid;
        type_  = vinfo->type;

        if (vinfo->ndim > 0)
        {
            const int ndim = vinfo->ndim;
            const int nblocks = vinfo->sum_nblocks;

            // get block-decomposition
            int err = adios_inq_var_blockinfo(file, vinfo);
            if (err)
            {
                WarningInFunction
                    << "Error reading blockinfo for dataset " << name_
                    << " from adios file: "
                    << adios_errmsg() << endl;
            }
            else
            {
                ADIOS_VARBLOCK *bp = vinfo->blockinfo;

                for (int blockI=0; blockI < nblocks; ++blockI)
                {
                    if (Foam::label(bp->process_id) == Pstream::myProcNo())
                    {
                        int count = 1;
                        for (int dimI=0; dimI < ndim; ++dimI)
                        {
                            count *= bp->count[dimI];
                        }
                        nElem_ += count;

#if 0
                        OStringStream os1;
                        OStringStream os2;

                        os1 << bp->count[0];
                        os2 << vinfo->dims[0];

                        for (int dimI=1; dimI < vinfo->ndim; ++dimI)
                        {
                            os1 << ',' << bp->count[dimI];
                            os2 << ',' << vinfo->dims[dimI];
                        }

                        if (vinfo->global)
                        {
                            Pout<< name_
                                << " [" << bp->process_id
                                << "] block[" << blockI << "]"
                                << " nElem: " << nElem_
                                << " dims: " << os1.str().c_str()
                                << " global: " << os2.str().c_str()
                                << " start: " << *(bp->start)
                                << " time-index:" << bp->time_index << endl;
                        }
                        else
                        {
                            Pout<< name_
                                << " [" << bp->process_id
                                << "] block[" << blockI << "]"
                                << " nElem: " << nElem_
                                << " dims: " << os1.str().c_str()
                                << " time-index:" << bp->time_index << endl;
                        }
#endif
                    }

                    ++bp;
                }
            }

            if (!nElem_)
            {
                WarningInFunction
                    << "No corresponding information for dataset " << name_
                    << " from adios file" << endl;
            }
        }
        else
        {
#if 0
            Pout<< " variable=" << name_ << " is scalar" << endl;
#endif
            nElem_ = 1; // scalar value
        }

        nBytes_ = nElem_ * adios_type_size(type_, NULL);

        if (verbose)
        {
            Pout<< " variable=" << name_
                << " nElem=" << nElem_
                << " nBytes=" << nBytes_
                << " type=" << adios_type_to_string(type_)
                << endl;
        }


        // free ADIOS_VARINFO and any ADIOS_VARBLOCK(s)
        adios_free_varinfo(vinfo);
    }

    return nBytes_ > 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adiosReader::VarInfo::VarInfo
(
    const char* varName
)
:
    name_(varName),
    varid_(0),
    type_(adios_unknown),
    nElem_(0),
    nBytes_(0)
{}


Foam::adiosReader::VarInfo::VarInfo
(
    const string& varName
)
:
    name_(varName),
    varid_(0),
    type_(adios_unknown),
    nElem_(0),
    nBytes_(0)
{}


Foam::adiosReader::VarInfo::VarInfo
(
    ADIOS_FILE *file,
    const char* varName,
    bool verbose
)
:
    name_(varName),
    varid_(0),
    type_(adios_unknown),
    nElem_(0),
    nBytes_(0)
{
    readInfo(file, verbose);
}


Foam::adiosReader::VarInfo::VarInfo
(
    ADIOS_FILE *file,
    const string& varName,
    bool verbose
)
:
    name_(varName),
    varid_(0),
    type_(adios_unknown),
    nElem_(0),
    nBytes_(0)
{
    readInfo(file, verbose);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::adiosReader::VarInfo::~VarInfo()
{}


// ************************************************************************* //
