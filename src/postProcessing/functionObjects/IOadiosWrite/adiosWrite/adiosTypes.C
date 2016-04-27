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

#include "adiosTypes.H"
#include "pTraits.H"

#include "adios.h"
#include "adios_read.h"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#if WM_LABEL_SIZE == 32

const enum ADIOS_DATATYPES
Foam::adiosTraits<Foam::label>::adiosType = adios_integer;

const size_t
Foam::adiosTraits<Foam::label>::adiosSize = 4;

#elif WM_LABEL_SIZE == 64

const enum ADIOS_DATATYPES
Foam::adiosTraits<Foam::label>::adiosType = adios_long;

const size_t
Foam::adiosTraits<Foam::label>::adiosSize = 8;

#else
# error "WM_LABEL_SIZE not defined in adiosTypes"
#endif

const int
Foam::adiosTraits<Foam::label>::nBits = WM_LABEL_SIZE;


#if defined(WM_SP)

const enum ADIOS_DATATYPES
Foam::adiosTraits<Foam::scalar>::adiosType = adios_real;

const size_t
Foam::adiosTraits<Foam::scalar>::adiosSize = 4;

const char* const
Foam::adiosTraits<Foam::scalar>::precisionName = "single";

#elif defined(WM_DP)

const enum ADIOS_DATATYPES
Foam::adiosTraits<Foam::scalar>::adiosType = adios_double;

const size_t
Foam::adiosTraits<Foam::scalar>::adiosSize = 8;

const char* const
Foam::adiosTraits<Foam::scalar>::precisionName = "double";

#else
#  error "WM_SP/WM_DP not defined in adiosTypes"
#endif


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::adiosTraits<Foam::label>::ok()
{
    if
    (
        sizeof(label) == adiosSize
     && adiosSize == adios_type_size(adiosType, NULL)
    )
    {
        return true;
    }
    else
    {
        FatalErrorInFunction
            << "Incorrect adios type/size for " << pTraits<label>::typeName
            << exit(FatalError);

        return false;
    }
}


bool Foam::adiosTraits<Foam::scalar>::ok()
{
    if
    (
        sizeof(scalar) == adiosSize
     && adiosSize == adios_type_size(adiosType, NULL)
    )
    {
        return true;
    }
    else
    {
        FatalErrorInFunction
            << "Incorrect adios type/size for " << pTraits<scalar>::typeName
            << exit(FatalError);

        return false;
    }
}


// ************************************************************************* //
