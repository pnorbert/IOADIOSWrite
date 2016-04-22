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
    Foam::IBufStream

Description
    A version of IStringStream that uses an externally managed buffer for
    its input.
    This allows the input buffer to be filled (and refilled) from
    various sources.

    Note that this stream will normally be used as a "one-shot" reader.
    Caution should be exercised that the referenced buffer remains valid
    (and ideally without any intermediate resizing) for the duration of
    the stream's use.

    An example of possible use:
    \code
        DynamicCharList buffer(4096);     // allocate some large buffer

        nread = something.read(buffer.data(),1024); // fill with content
        buffer.setSize(nread);              // content size

        // construct dictionary, or do something else
        IBufStream is(buffer)
        dictionary dict1(is);

        // sometime later
        nread = something.read(buffer.data(),2048); // fill with content
        buffer.setSize(nread);              // content size

        // without intermediate variable
        dictionary dict2(IBufStream(buffer)());
    \endcode
\*---------------------------------------------------------------------------*/

#ifndef IBufStream_H
#define IBufStream_H

#include "IStringStream.H"
#include "DynamicCharList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class IBufStream Declaration
\*---------------------------------------------------------------------------*/

class IBufStream
:
    public IStringStream
{
private:
    // Private data

        //- Reference to the underlying character buffer
        UList<char>& buf_;

public:

    // Constructors

        //- Construct using specified buffer for the input
        IBufStream
        (
            const UList<char>& buffer,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            IStringStream(format, version),
            buf_(const_cast<UList<char>&>(buffer))
        {
            // Associate character buffer with stream
            stdStream().rdbuf()->pubsetbuf(buf_.data(), buf_.size());
        }


    // Member operators

        //- Return a non-const reference to const Istream
        //  Needed for read-constructors where the stream argument is temporary:
        //  e.g. thing thisThing(IFstream("thingFileName")());
        Istream& operator()() const
        {
            return const_cast<Istream&>(static_cast<const Istream&>(*this));
        }

};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
