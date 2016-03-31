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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::adiosWrite::fieldDefine(label regionID)
{
    const regionInfo& rInfo = regions_[regionID];

    char tmpstr[128];
    Info<< "  adiosWrite::fieldDefine: region " << regionID << " "
        << rInfo.name_ << ": " << endl;

    sprintf(tmpstr, "region%d", regionID);
    adios_define_attribute
    (
        groupID_,
        "name",
        tmpstr,
        adios_string,
        regions_[regionID].name_.c_str(),
        NULL
    );

    const fvMesh& mesh = time_.lookupObject<fvMesh>(regions_[regionID].name_);
    fieldDefine<volScalarField>(mesh, rInfo.scalarFields_, regionID);
    fieldDefine<volVectorField>(mesh, rInfo.vectorFields_, regionID);
    fieldDefine<surfaceScalarField>(mesh, rInfo.surfaceScalarFields_, regionID);
}


void Foam::adiosWrite::fieldWrite(label regionID)
{
    const regionInfo& rInfo = regions_[regionID];

    Info<< "  adiosWrite::fieldWrite: region " << regionID << " "
        << rInfo.name_ << ": " << endl;

    const fvMesh& mesh = time_.lookupObject<fvMesh>(regions_[regionID].name_);
    fieldWrite<volScalarField>(mesh, rInfo.scalarFields_, regionID);
    fieldWrite<volVectorField>(mesh, rInfo.vectorFields_, regionID);
    fieldWrite<surfaceScalarField>(mesh, rInfo.surfaceScalarFields_, regionID);
}


// ************************************************************************* //
