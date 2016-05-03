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
    Similar to IStringStream but using an externally managed buffer for its
    input. This allows the input buffer to be filled (and refilled) from
    various sources.

    Note that this stream will normally be used as a "one-shot" reader.
    Caution must be exercised that the referenced buffer remains valid and
    without any intermediate resizing for the duration of the stream's use.

    An example of possible use:
    \code
        DynamicList<char> buffer(4096);   // allocate some large buffer

        nread = something.read(buffer.data(),1024); // fill with content
        buffer.setSize(nread);            // content size

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

#include "OBufStream.H"

#include "ISstream.H"
#include "UList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class IBufStreamAllocator Declaration
\*---------------------------------------------------------------------------*/

//- An allocator for std::istream with external buffers
class IBufStreamAllocator
{
private:

    // Private data

        //- Reference to the underlying buffer
        memorybuf::input buf_;

protected:

    // Protected data

        //- The stream pointer
        std::istream* strPtr_;

public:

    // Constructors

        //- Construct with buffer and number of bytes
        IBufStreamAllocator(const char *buffer, size_t nbytes)
        :
            buf_(buffer, nbytes),
            strPtr_(new std::istream(&buf_))
        {}


        //- Construct with buffer and number of bytes
        IBufStreamAllocator(const void *buffer, size_t nbytes)
        :
            buf_(buffer, nbytes),
            strPtr_(new std::istream(&buf_))
        {}


        //- Construct using specified buffer and its inherent storage size
        template<class T>
        IBufStreamAllocator(const UList<T>& buffer)
        :
            buf_(buffer),
            strPtr_(new std::istream(&buf_))
        {}


    //- Destructor
    ~IBufStreamAllocator()
    {
        delete strPtr_;
    }


    // Member functions

        //- Return the current get position in the buffer
        std::streampos pos()
        {
            return strPtr_->tellg();
        }


        //- Move to buffer start, clear errors
        void rewind()
        {
            strPtr_->rdbuf()->pubseekpos(0);
            strPtr_->clear(); // for safety, clear any old errors
        }

};


/*---------------------------------------------------------------------------*\
                         Class IBufStream Declaration
\*---------------------------------------------------------------------------*/

class IBufStream
:
    private IBufStreamAllocator,
    public ISstream
{
public:

    // Constructors

        //- Construct using specified buffer and number of bytes
        IBufStream
        (
            const char* buffer,
            size_t nbytes,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            IBufStreamAllocator(buffer, nbytes),
            ISstream
            (
                *strPtr_,
                "IBufStream.source",
                format,
                version
            )
        {}


        //- Construct using data area from buffer and number of bytes
        template<class T>
        IBufStream
        (
            const UList<T>& buffer,
            size_t nbytes,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            IBufStreamAllocator(buffer.cdata(), nbytes),
            ISstream
            (
                *strPtr_,
                "IBufStream.source",
                format,
                version
            )
        {}


        //- Construct using specified buffer and its inherent storage size
        //  Uses addressed size, thus no special treatment for a DynamicList
        template<class T>
        IBufStream
        (
            const UList<T>& buffer,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            IBufStreamAllocator(buffer),
            ISstream
            (
                *strPtr_,
                "IBufStream.source",
                format,
                version
            )
        {}


    //- Destructor
    ~IBufStream()
    {}


    // Member functions

        //- Return the current input position in the buffer
        using IBufStreamAllocator::pos;


        //- Rewind the stream, clearing any old errors
        virtual Istream& rewind()
        {
            IBufStreamAllocator::rewind();
            setGood();  // resynchronize with internal state

            return *this;
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
