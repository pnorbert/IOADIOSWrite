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
    Foam::ParticleBinaryBlob

Description
    Access routines for addressing parcel/particle binary content.

\*---------------------------------------------------------------------------*/

#ifndef ParticleBinaryBlob_H
#define ParticleBinaryBlob_H

#include "DynamicList.H"
#include "ParticleBinaryBlobFragment.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   // forward declaration
   class ParticleBinaryBlob;


/*---------------------------------------------------------------------------*\
                   Class ParticleBinaryFragment Declaration
\*---------------------------------------------------------------------------*/
class ParticleBinaryBlob
:
    private SLList<ParticleBinaryBlobFragment>
{
public:

    typedef SLList<ParticleBinaryBlobFragment> Container;
    typedef ParticleBinaryBlobFragment Fragment;

private:

    // Private data

        wordList types_;
        wordList names_;


    // Private Member Functions

        //- Clear internal lists
        void clear();

        //- Set fragment types and sizes, add names later
        void setTypes
        (
            const UList<word>& inputTypes,
            const bool raw
        );

        //- Set fragment types and sizes, add names later
        void setTypes
        (
            Istream& is,
            const bool raw
        );

        //- Set fragment types and sizes, add names later
        void setTypes
        (
            const string& inputTypes,
            const bool raw
        );

        //- Set the names for each fragment
        void setNames(const UList<word>& inputNames);

        //- Set the names for each fragment
        void setNames(Istream& is);

        //- Set the names for each fragment
        void setNames(const string& inputNames);

        //- Fill in lookup tables
        void makeSummary();

public:

    // Constructors

        //- Construct null
        ParticleBinaryBlob();

        //- Construct from base components (types, names, offsets, widths)
        ParticleBinaryBlob
        (
            const UList<word>& inputTypes,
            const UList<word>& inputNames,
            const UList<int>& offsetList,
            const UList<int>& bytesList
        );

        //- Construct from type names and variable names.
        //  Optionally specify that the expected content is 'raw'
        //  (ie, without any surrounding brackets)
        ParticleBinaryBlob
        (
            const UList<word>& inputTypes,
            const UList<word>& inputNames,
            const bool raw = false
        );

        //- Construct from type names and variable names.
        //  Optionally specify that the expected content is 'raw'
        //  (ie, without any surrounding brackets)
        ParticleBinaryBlob
        (
            const string& inputTypes,
            const string& inputNames,
            const bool raw = false
        );


        //- Return the binary size in number of characters of the UList
        //  if the element is a primitive type
        //  i.e. contiguous<T>() == true.
        //  Note that is of type streamsize since used in stream ops
        inline size_t byteSize() const
        {
            if (empty())
            {
                return 0;
            }

            // The overall size must include a closing ')' character.
            // If the first fragment has a non-zero offset it has an
            // opening '(' character
            return last().end() + (first().offset() ? 1 : 0);
        }


        //- A summary of the fragment types
        const wordList& types() const;

        //- A summary of the fragment names
        const wordList& names() const;

        //- A summary of the fragment offsets
        List<int> offsets() const;

        //- A summary of the fragment widths
        List<int> sizes() const;


        using Container::size;
        using Container::begin;
        using Container::end;
        using Container::cbegin;
        using Container::cend;
        using Container::first;
        using Container::last;


    // IOstream operators

        friend Ostream& operator<<
        (
            Ostream& os,
            const ParticleBinaryBlob& lst
        );
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
