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

Class
    Foam::ParticleBinaryBlobFragment

Description
    Access for fragments of parcel/particle binary content.

\*---------------------------------------------------------------------------*/

#ifndef ParticleBinaryBlobFragment_H
#define ParticleBinaryBlobFragment_H

#include "wordList.H"
#include "label.H"
#include "scalar.H"
#include "vector.H"

#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   // forward declaration
   class ParticleBinaryBlob;

/*---------------------------------------------------------------------------*\
                 Class ParticleBinaryBlobFragment Declaration
\*---------------------------------------------------------------------------*/
class ParticleBinaryBlobFragment
{
   friend ParticleBinaryBlob;

    // Private data

        //- Offset in bytes from beginning of binary blob
        const size_t offset_;

        //- Width of primitive type
        size_t width_;

        //- Number of primitive types
        size_t count_;

        //- Number of components for the primitive type
        int nComponents;

        //- Name of enclosing container, if any (future use)
        word container_;

        //- Name of the primitive (base) type
        word type_;

        //- Name of the associated variable
        word name_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const ParticleBinaryBlobFragment&);


        //- Width and nComponents for supported primitives
        template<class PrimitiveType>
        bool setPrimitiveType()
        {
            width_      = sizeof(pTraits<PrimitiveType>);
            nComponents = pTraits<PrimitiveType>::nComponents;
            return true;
        }


        //- Width for supported primitives
        static size_t width(const word& typeTok);

        //- Width and nComponents for supported primitives
        bool set(const word& typeTok);


public:

    //- Construct for specified type and offset
    ParticleBinaryBlobFragment
    (
        const word& type,
        size_t offset
    );


    //- The name associated with this byte sequence
    inline const word& name() const
    {
        return name_;
    }


    //- The primitive type for this byte sequence
    inline const word& type() const
    {
        return type_;
    }


    //- The beginning of this byte sequence
    inline size_t offset() const
    {
        return offset_;
    }


    //- The size of an individual element
    inline size_t width() const
    {
        return width_;
    }


    //- The number of individual elements
    inline size_t count() const
    {
        return count_;
    }


    //- Total number of bytes
    inline size_t sizeOf() const
    {
        return width_ * count_;
    }


    //- The end of this byte sequence
    inline size_t end() const
    {
        return offset_ + sizeOf();
    }



    // TODO
    // restrict to the respective type

    inline label getLabel(const char* buffer) const
    {
        return *(reinterpret_cast<const label*>(buffer + offset_));
    }


    inline scalar getScalar(const char* buffer) const
    {
        return *(reinterpret_cast<const scalar*>(buffer + offset_));
    }


    inline vector getVector(const char* buffer) const
    {
        return *(reinterpret_cast<const vector*>(buffer + offset_));
    }


    inline label getLabel(const UList<char>& buffer) const
    {
        return getLabel(buffer.cdata());
    }


    inline scalar getScalar(const UList<char>& buffer) const
    {
        return getScalar(buffer.cdata());
    }


    inline vector getVector(const UList<char>& buffer) const
    {
        return getVector(buffer.cdata());
    }


// IOstream operators

    friend Ostream& operator<<
    (
        Ostream& os,
        const ParticleBinaryBlobFragment&
    );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
