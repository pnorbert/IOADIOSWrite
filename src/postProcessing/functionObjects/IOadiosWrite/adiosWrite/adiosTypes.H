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

Typedef
    Foam::Scalar

Description
    Single floating point number (float or double)

SourceFiles
    Scalar.C

\*---------------------------------------------------------------------------*/

#include "adiosTraits.H"

#include "label.H"
#include "scalar.H"
#include "vector.H"

#include "adios_types.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- template specialization for adiosTraits<label>
template<>
class adiosTraits<label>
{
    label p_;

public:

    //- Component type
    typedef label cmptType;

    //- ADIOS type
    static const enum ADIOS_DATATYPES adiosType;

    //- Size of ADIOS type
    static const size_t adiosSize;

    //- Number of bits
    static const int nBits;

    //- Verify ADIOS type-size
    static bool ok();


    // Constructors

        //- Construct from primitive
        explicit adiosTraits(const label&);


    // Member Functions

        //- Access to the primitive value
        operator label() const
        {
            return p_;
        }

        //- Access to the primitive value
        operator label&()
        {
            return p_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- template specialization for adiosTraits<scalar>
template<>
class adiosTraits<scalar>
{
    scalar p_;

public:

    //- Component type
    typedef scalar cmptType;

    //- ADIOS type
    static const enum ADIOS_DATATYPES adiosType;

    //- Size of ADIOS type
    static const size_t adiosSize;

    //- Name of the precision (single|double)
    static const char* const precisionName;

    //- Verify ADIOS type-size
    static bool ok();


    // Constructors

        //- Construct from primitive
        explicit adiosTraits(const scalar&);


    // Member Functions

        //- Access to the primitive value
        operator scalar() const
        {
            return p_;
        }

        //- Access to the primitive value
        operator scalar&()
        {
            return p_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
