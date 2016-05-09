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

#include "adios.h"
#include "adios_read.h"


// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //

template<class T>
size_t Foam::adiosReader::getBuffered
(
    const string& name,
    DynamicList<T>& buffer
) const
{
    size_t nbytes = 0;

    if (hasVariable(name))
    {
        const VarInfo& vinfo = variables[name];
        nbytes = vinfo.sizeOf();

        const size_t sizeT = sizeof(T);
        const size_t target = (nbytes / sizeT) + (nbytes % sizeT);

        buffer.reserve(nbytes / sizeT + nbytes % sizeT);

        if (getVariable(name, buffer.data()))
        {
            buffer.setSize(target);
        }
        else
        {
            buffer.setSize(0);
            nbytes = 0;
        }
    }

    if (!nbytes)
    {
        FatalErrorInFunction
            << "missing adios variable: " << name
            << exit(FatalIOError);
    }

    return nbytes;
}



template<class StringType>
bool Foam::adiosReader::readStringListAttributeIfPresent
(
    const string& name,
    List<StringType>& lst
) const
{
    enum ADIOS_DATATYPES type = adios_unknown;
    int  size = 0;
    void *data = 0;

    bool ok =
    (
        hasAttribute(name)
     && 0 ==
        adios_get_attr_byid
        (
            file,
            attributes[name],
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
                << "string-list attribute has wrong type: " << name << nl
                << "  expecting char*[], found "
                << adios_type_to_string(type)
                << exit(FatalIOError);
        }
    }

    return ok;
}


// ************************************************************************* //
