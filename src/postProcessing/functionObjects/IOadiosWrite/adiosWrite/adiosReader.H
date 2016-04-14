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

Class
    Foam::adiosReader

Description

SourceFiles
    adiosReader.C

\*---------------------------------------------------------------------------*/

#ifndef adiosReader_H
#define adiosReader_H

#include "fileName.H"
#include "IStringStreamBuf.H"

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

        ADIOS_FILE *file;
        ADIOS_SELECTION *selection;
        size_t maxLen;

        helper();

        ~helper();


        bool open(const fileName& bpFile, MPI_Comm comm);

        void select(ADIOS_SELECTION *sel);

        bool getDataSet(const fileName& varName, void* data);

        bool getDataSet(const fileName& varName, IStringStreamBuf& is);

        void scan(bool verbose=false);

        void close();

    };

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
