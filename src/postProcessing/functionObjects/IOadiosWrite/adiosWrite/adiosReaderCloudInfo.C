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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::adiosReader::cloudInfo::cloudInfo()
:
    fileName(),
    className_(),
    nBytes_(0),
    nTotal_(0),
    nParticle_(0),
    width_(0)
{}


Foam::adiosReader::cloudInfo::cloudInfo(const fileName& varName)
:
    fileName(varName),
    className_(),
    nBytes_(0),
    nTotal_(0),
    nParticle_(0),
    width_(0)
{}


Foam::adiosReader::cloudInfo::cloudInfo
(
    const adiosReader::VarInfo& varinfo,
    const adiosReader& reader
)
:
    fileName(varinfo.name()),
    className_(),
    nBytes_(varinfo.sizeOf()),
    nTotal_(0),
    nParticle_(0),
    width_(0)
{
    read(reader);
}


bool Foam::adiosReader::cloudInfo::read
(
    const adiosReader& reader
)
{
    const fileName& varName = fullName();

    nParticle_ = reader.getIntVariable(varName/"nParticle");

    className_ = reader.getStringAttribute(varName/"class");
    nTotal_    = reader.getIntAttribute(varName/"nParticle");
    width_     = reader.getIntAttribute(varName/"size");

    blob_ = ParticleBinaryBlob
    (
        reader.getStringListAttribute<word>(varName/"types"),
        reader.getStringListAttribute<word>(varName/"names"),
        false, // raw = uncooked
        reader.getIntListAttribute(varName/"offset"),
        reader.getIntListAttribute(varName/"byte-size")
    );

    return nBytes_ > 0;
}


// ************************************************************************* //
