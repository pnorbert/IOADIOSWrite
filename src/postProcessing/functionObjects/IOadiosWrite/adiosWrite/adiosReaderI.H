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

#include "fileName.H"
#include "IBufStream.H"
#include "DynamicCharList.H"
#include "HashTable.H"

#include "Ostream.H"

#include "adios.h"
#include "adios_read.h"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline const Foam::fileName& Foam::adiosReader::VarInfo::name() const
{
    return name_;
}


inline enum ADIOS_DATATYPES Foam::adiosReader::VarInfo::dataType() const
{
    return type_;
}


inline Foam::label Foam::adiosReader::VarInfo::nElem() const
{
    return nElem_;
}


size_t Foam::adiosReader::VarInfo::sizeOf() const
{
    return nBytes_;
}


inline size_t Foam::adiosReader::sizeOf() const
{
    return maxLen;
}


inline bool Foam::adiosReader::hasAttribute(const string& name) const
{
    return attributes.found(name);
}


inline bool Foam::adiosReader::hasVariable(const string& name) const
{
    return variables.found(name);
}


inline Foam::label
Foam::adiosReader::getIntAttribute(const string& name) const
{
    label value;

    if (!readIntAttributeIfPresent(name, value))
    {
        FatalErrorInFunction
            << "integer attribute missing: " << name
            << exit(FatalIOError);
    }

    return value;
}


inline Foam::List<Foam::label>
Foam::adiosReader::getIntListAttribute(const string& name) const
{
    List<label> value;

    if (!readIntListAttributeIfPresent(name, value))
    {
        FatalErrorInFunction
            << "int-list attribute missing: " << name
            << exit(FatalIOError);
    }

    return value;
}


inline double
Foam::adiosReader::getScalarAttribute(const string& name) const
{
    double value;

    if (!readScalarAttributeIfPresent(name, value))
    {
        FatalErrorInFunction
            << "double attribute missing: " << name
            << exit(FatalIOError);
    }

    return value;
}


inline Foam::List<double>
Foam::adiosReader::getScalarListAttribute(const string& name) const
{
    List<double> value;

    if (!readScalarListAttributeIfPresent(name, value))
    {
        FatalErrorInFunction
            << "double-list attribute missing: " << name
            << exit(FatalIOError);
    }

    return value;
}


inline Foam::string
Foam::adiosReader::getStringAttribute(const string& name) const
{
    string value;

    if (!readStringAttributeIfPresent(name, value))
    {
        FatalErrorInFunction
            << "string attribute missing: " << name
            << exit(FatalIOError);
    }

    return value;
}


inline Foam::label
Foam::adiosReader::getIntVariable(const string& name) const
{
    label value;

    if (!readIntVariableIfPresent(name, value))
    {
        FatalErrorInFunction
            << "integer variable missing: " << name
            << exit(FatalIOError);
    }

    return value;
}


inline Foam::scalar
Foam::adiosReader::getScalarVariable(const string& name) const
{
    scalar value;

    if (!readScalarVariableIfPresent(name, value))
    {
        FatalErrorInFunction
            << "scalar value missing: " << name
            << exit(FatalIOError);
    }

    return value;
}


template<class StringType>
inline Foam::List<StringType> Foam::adiosReader::getStringListAttribute
(
    const string& name
) const
{
    List<StringType> value;

    if (!readStringListAttributeIfPresent(name, value))
    {
        FatalErrorInFunction
            << "string-list attribute missing: " << name
            << exit(FatalIOError);
    }

    return value;
}


// ************************************************************************* //
