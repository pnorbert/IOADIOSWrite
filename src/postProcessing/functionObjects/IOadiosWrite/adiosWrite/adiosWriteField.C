/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |               2015 Norbert Podhorszki
                            |               2016 OpenCFD Ltd.
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
#include "OStringStream.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

size_t Foam::adiosWrite::fieldDefine(label regionID)
{
    const regionInfo& rInfo = regions_[regionID];

    Info<< "  adiosWrite::fieldDefine: region " << regionID << "="
        << rInfo.name_ << endl;

    adios_define_attribute
    (
        groupID_,
        "name",
        ("region" + Foam::name(regionID)).c_str(),
        adios_string,
        regions_[regionID].name_.c_str(),
        NULL
    );

    const fvMesh& mesh = time_.lookupObject<fvMesh>(regions_[regionID].name_);

    size_t bufLen = 0;
    size_t maxLen = 0;

    bufLen = fieldDefine<volScalarField>(mesh, rInfo.scalarFields_, regionID);
    maxLen = Foam::max(maxLen, bufLen);

    bufLen = fieldDefine<volVectorField>(mesh, rInfo.vectorFields_, regionID);
    maxLen = Foam::max(maxLen, bufLen);

    bufLen = fieldDefine<surfaceScalarField>(mesh, rInfo.surfaceScalarFields_, regionID);
    maxLen = Foam::max(maxLen, bufLen);

    return maxLen;
}


void Foam::adiosWrite::fieldWrite(label regionID)
{
    const regionInfo& rInfo = regions_[regionID];
    const fvMesh& mesh = time_.lookupObject<fvMesh>(rInfo.name_);

    Info<< "  adiosWrite::fieldWrite: region " << regionID << "="
        << rInfo.name_ << endl;

    fieldWrite<volScalarField>(mesh, rInfo.scalarFields_, regionID);
    fieldWrite<volVectorField>(mesh, rInfo.vectorFields_, regionID);
    fieldWrite<surfaceScalarField>(mesh, rInfo.surfaceScalarFields_, regionID);
}


// ************************************************************************* //
