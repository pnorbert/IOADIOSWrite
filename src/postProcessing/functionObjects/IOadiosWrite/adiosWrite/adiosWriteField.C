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

#include "adiosWrite.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

size_t Foam::adiosWrite::fieldDefine(const regionInfo& rInfo)
{
    Info<< "  adiosWrite::fieldDefine: " << rInfo.info() << endl;

    const fvMesh& mesh = time_.lookupObject<fvMesh>(rInfo.name_);

    size_t bufLen = 0;
    size_t maxLen = 0;

    bufLen = fieldDefine<volScalarField>(mesh, rInfo, rInfo.scalarFields_);
    maxLen = Foam::max(maxLen, bufLen);

    bufLen = fieldDefine<volVectorField>(mesh, rInfo, rInfo.vectorFields_);
    maxLen = Foam::max(maxLen, bufLen);

    bufLen = fieldDefine<surfaceScalarField>(mesh, rInfo, rInfo.surfaceScalarFields_);
    maxLen = Foam::max(maxLen, bufLen);

    return maxLen;
}


void Foam::adiosWrite::fieldWrite(const regionInfo& rInfo)
{
    Info<< "  adiosWrite::fieldWrite: " << rInfo.info() << endl;

    const fvMesh& mesh = time_.lookupObject<fvMesh>(rInfo.name_);

    fieldWrite<volScalarField>(mesh, rInfo, rInfo.scalarFields_);
    fieldWrite<volVectorField>(mesh, rInfo, rInfo.vectorFields_);
    fieldWrite<surfaceScalarField>(mesh, rInfo, rInfo.surfaceScalarFields_);
}


// ************************************************************************* //
