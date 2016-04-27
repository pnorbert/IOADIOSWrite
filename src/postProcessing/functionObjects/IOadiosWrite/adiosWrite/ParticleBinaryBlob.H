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

    typedef SLList<ParticleBinaryBlobFragment> Container;

public:

    typedef ParticleBinaryBlobFragment Fragment;

private:

    // Private data

        wordList types_;
        wordList names_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        ParticleBinaryBlob(const ParticleBinaryBlob&);

        //- Disallow default bitwise assignment
        void operator=(const ParticleBinaryBlob&);


        //- Clear internal lists
        void clear();

        //- Set fragment types and sizes, add names later
        void setTypes(const UList<word>& inputTypes);

        //- Set fragment types and sizes, add names later
        void setTypes(Istream& is);

        //- Set fragment types and sizes, add names later
        void setTypes(const string& inputTypes);

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

        //- Construct from type names and variable names
        //  Optionally verify the offsets and widths
        ParticleBinaryBlob
        (
            const UList<word>& inputTypes,
            const UList<word>& inputNames,
            const UList<int>& offsetList = UList<int>::null(),
            const UList<int>& bytesList  = UList<int>::null()
        );


        //- Construct from type names and variable names
        ParticleBinaryBlob
        (
            const string& inputTypes,
            const string& inputNames
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
            else
            {
                return last().end();
            }
        }


        const wordList& types() const;

        const wordList& names() const;

        List<int> offsets() const;

        List<int> byteSizes() const;


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
