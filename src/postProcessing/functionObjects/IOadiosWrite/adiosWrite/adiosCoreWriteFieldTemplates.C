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

#include "adiosCoreWrite.H"

#include "OBufStream.H"
#include "IOstream.H"
#include "IOstreams.H"
#include "Ostream.H"


// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //

template<class FieldType>
int64_t Foam::adiosCoreWrite::defineInternalField
(
    const FieldType& field,
    const fileName& basePath
)
{
    const fileName varName = basePath/field.name();

    // could also write field.dimensions() as an attribute (as required)
    defineAttribute("class", varName, field.type());

    OutputCounter counter(adiosCore::strFormat);
    counter << field;

    return defineStreamVariable(varName, counter.size());
}


template<class FieldType>
int64_t Foam::adiosCoreWrite::defineField
(
    const FieldType& field,
    const fileName& basePath
)
{
    const fileName varName = basePath/field.name();

    const typename FieldType::GeometricBoundaryField& bfield =
        field.boundaryField();

    // independent of how we store fields,
    // a quick lookup of field patch types may prove useful
    stringList pTypes(bfield.size());
    forAll(bfield, patchI)
    {
        pTypes[patchI] = bfield[patchI].type();
    }
    defineListAttribute("patch-types", varName, pTypes);

    return defineInternalField(field, basePath);
}


template<class FieldType>
void Foam::adiosCoreWrite::writeField
(
    const FieldType& field,
    const fileName& basePath
)
{
    // use iobuffer_ to avoid reallocations
    // Needs rework?
    iobuffer_.reserve(adiosCoreWrite::maxSize());

    OutputBufStreamer os(iobuffer_, adiosCore::strFormat);
    os << field;

    // Do the actual write (as stream)
    writeVariable(basePath/field.name(), iobuffer_);
}


// ************************************************************************* //
