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

#include "ParticleBinaryBlobFragment.H"

#include "vector.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

size_t Foam::ParticleBinaryBlobFragment::width(const word& typeTok)
{
    if (typeTok == pTraits<label>::typeName)
    {
        return sizeof(label);
    }
    if (typeTok == pTraits<scalar>::typeName)
    {
        return sizeof(scalar);
    }
    if (typeTok == pTraits<vector>::typeName)
    {
        return sizeof(vector);
    }

    return 0;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::ParticleBinaryBlobFragment::set(const word& typeTok)
{
    if (typeTok == pTraits<label>::typeName)
    {
        return setPrimitiveType<label>();
    }
    if (typeTok == pTraits<scalar>::typeName)
    {
        return setPrimitiveType<scalar>();
    }
    if (typeTok == pTraits<vector>::typeName)
    {
        return setPrimitiveType<vector>();
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ParticleBinaryBlobFragment::ParticleBinaryBlobFragment
(
    const word& type,
    size_t offset
)
:
    offset_(offset),
    width_(0),
    count_(0),
    nComponents(1),
    container_(),
    type_(type),
    name_("__unnamed__")
{
    set(type_);

    width_ = width(type_);

    // Pair<..>
    if (width_ > 0)
    {
        count_ = 1;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ParticleBinaryBlobFragment& frag
)
{
    // Write contents

    os.writeKeyword("type")   << frag.type_     << token::END_STATEMENT << endl;
    os.writeKeyword("name")   << frag.name_     << token::END_STATEMENT << endl;
    os.writeKeyword("offset") << frag.offset()  << token::END_STATEMENT << endl;
    os.writeKeyword("width")  << frag.width()   << token::END_STATEMENT << endl;
    os.writeKeyword("count")  << frag.count()   << token::END_STATEMENT << endl;
    os.writeKeyword("end")    << frag.end()     << token::END_STATEMENT << endl;

    os.check
    (
        "Ostream& operator<<(Ostream&, const ParticleBinaryBlob::Fragment&)"
    );

    return os;
}


// ************************************************************************* //
