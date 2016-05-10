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

#include "adiosTime.H"
#include "adiosReader.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::adiosTime::attrType,
        4
    >::names[] =
    {
        "/time/index",
        "/time/value",
        "/time/deltaT",
        "/time/deltaT0"
    };
}


const Foam::NamedEnum<Foam::adiosTime::attrType, 4>
Foam::adiosTime::attr;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adiosTime::adiosTime()
:
    index_(-1),
    value_(VGREAT),
    deltaT_(0),
    deltaT0_(0)
{}


Foam::adiosTime::adiosTime(const TimeState& t)
:
    index_(t.timeIndex()),
    value_(t.timeOutputValue()),
    deltaT_(t.deltaT().value()),
    deltaT0_(t.deltaT0().value())
{}


Foam::adiosTime::adiosTime(const adiosReader& reader)
:
    index_(reader.getIntVariable(attr[INDEX])),
    value_(reader.getScalarVariable(attr[VALUE])),
    deltaT_(reader.getScalarVariable(attr[DT])),
    deltaT0_(reader.getScalarVariable(attr[DT0]))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::adiosTime::valid() const
{
    return index_ > 0;
}


// ************************************************************************* //
