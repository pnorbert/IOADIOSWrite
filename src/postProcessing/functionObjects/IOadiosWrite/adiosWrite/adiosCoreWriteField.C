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
#include "nullObject.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::adiosCoreWrite::supportedFieldType(const word& fieldType)
{
    if (isNull(fieldType) || fieldType.empty())
    {
        return false;
    }

    return
    (
        fieldType == volScalarField::typeName
     || fieldType == volVectorField::typeName
     || fieldType == surfaceScalarField::typeName
     || fieldType == volSphericalTensorField::typeName
     || fieldType == volSymmTensorField::typeName
     || fieldType == volTensorField::typeName

        // internal fields
     || fieldType == volScalarField::DimensionedInternalField::typeName
     || fieldType == volVectorField::DimensionedInternalField::typeName
    );
}


// ************************************************************************* //
