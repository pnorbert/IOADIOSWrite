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

#include "error.H"
#include "CStringList.H"

#include "Ostream.H"

#include "adios.h"
#include "adios_read.h"


// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //

template<class StringType>
bool Foam::adiosReader::helper::readStringListAttributeIfPresent
(
    const string& attrName,
    List<StringType>& lst
) const
{
    enum ADIOS_DATATYPES type = adios_unknown;
    int  size = 0;
    void *data = 0;

    bool ok =
    (
        attributes.found(attrName)
     && 0 ==
        adios_get_attr_byid
        (
            file,
            attributes[attrName],
            &type,
            &size,
            &data
        )
    );

    if (ok)
    {
        if (type == adios_string_array)
        {
            lst = CStringList::asList<StringType>
            (
                (size / sizeof(char*)),         // argc
                reinterpret_cast<char**>(data)  // argv
            );
        }
        else
        {
            ok = false;
        }

        if (data)
        {
            free(data);
        }

        if (!ok)
        {
            FatalErrorInFunction
                << "string-list attribute has wrong type: " << attrName << nl
                << "  expecting char*[], found "
                << adios_type_to_string(type)
                << exit(FatalIOError);
        }
    }

    return ok;
}


template<class StringType>
Foam::List<StringType> Foam::adiosReader::helper::getStringListAttribute
(
    const string& attrName
) const
{
    List<StringType> value;

    if (!readStringListAttributeIfPresent(attrName, value))
    {
        FatalErrorInFunction
            << "string-list attribute missing: " << attrName
            << exit(FatalIOError);
    }

    return value;
}


// ************************************************************************* //
