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
    Foam::OBufStream

Description
    A version of OStringStream that uses an externally managed buffer for
    its output.

    This allows the output buffer to be reused and can make it easier when
    writing out data.  It is the user's responsibility to ensure proper
    synchronization in the sizes. Provided that the external buffer is large
    enough that overflow does not occur, the following usage pattern
    should work.

    \code
        DynamicCharList buffer(4096);     // allocate some large buffer

        {
            OBufStream os(buffer);
            os << "content1" << " and more content";
            buffer.setSize(os.size());      // synchronize sizes
        }

        something.write(buffer, buffer.size());
    \endcode

    Although the OBufStream is quite lightweight, there may be cases
    where it is preferable to reuse the stream as well.
    \code
        DynamicCharList buffer(4096);     // allocate some large buffer

        OBufStream os(buffer);
        os << "content1" << " and more content";
        buffer.setSize(os.size());      // synchronize sizes

        something.write(buffer, buffer.size());

        os.rewind();
        os << "content2";
        buffer.setSize(os.size());      // synchronize sizes

        something.write(buffer, buffer.size());

        // or simply using the output size directly (without sync)
        os.rewind();
        os << "content3";

        something.write(buffer, os.size());
    \endcode

\*---------------------------------------------------------------------------*/

#ifndef OBufStream_H
#define OBufStream_H

#include "OStringStream.H"
#include "DynamicCharList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class OBufStreamAllocator Declaration
\*---------------------------------------------------------------------------*/

//- An allocator for std::ostream with external buffers
class OBufStreamAllocator
{
private:

    // Private data

        //- Reference to the underlying character buffer
        UList<char>& buf_;

protected:

    // Protected data

        //- The stream pointer
        ostream* strPtr_;

public:

    // Constructors

        //- Construct with buffer
        OBufStreamAllocator(const UList<char>& buffer)
        :
            buf_(const_cast<UList<char>&>(buffer)),
            strPtr_(new std::ostringstream())
        {
            strPtr_->rdbuf()->pubsetbuf(buf_.data(), buf_.size());
        }


        //- Construct with buffer and size
        OBufStreamAllocator(const UList<char>& buffer, label listSize)
        :
            buf_(const_cast<UList<char>&>(buffer)),
            strPtr_(new std::ostringstream())
        {
            strPtr_->rdbuf()->pubsetbuf(buf_.data(), listSize);
        }


    //- Destructor
    ~OBufStreamAllocator()
    {
        delete strPtr_;
    }


    //- Return the current output position in the buffer
    //  Providing that the external buffer is large enough that overflow does
    //  not occur, this should correspond to the number of bytes output.
    std::streampos size()
    {
        return strPtr_->tellp();
    }


    //- Rewind the output buffer
    void rewind()
    {
        strPtr_->rdbuf()->pubseekpos(0);
    }

};


/*---------------------------------------------------------------------------*\
                         Class OBufStream Declaration
\*---------------------------------------------------------------------------*/

//- A std::ostream with an attached external buffer
//
//  Providing that the external buffer is large enough that no overflow occurs,
//  it should be possible to determine the number of bytes output by querying
//  the output position in the buffer.
class OBufStream
:
    private OBufStreamAllocator,
    public OSstream
{
public:

    // Constructors

        //- Construct using specified buffer for the output
        OBufStream
        (
            const UList<char>& buffer,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            OBufStreamAllocator(buffer),
            OSstream
            (
                *strPtr_,
                "OStringStream.sinkFile", // keep same name as string-stream
                format,
                version
            )
        {}


        //- Construct using specified buffer for the output
        OBufStream
        (
            const DynamicCharList& buffer,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            OBufStreamAllocator(buffer, buffer.capacity()),
            OSstream
            (
                *strPtr_,
                "OStringStream.sinkFile", // keep same name as string-stream
                format,
                version
            )
        {}


    //- Destructor
    ~OBufStream()
    {}


    // Member functions

    //- Return the current output position in the buffer
    //  Providing that the external buffer is large enough that overflow does
    //  not occur, this should correspond to the number of bytes output.
    using OBufStreamAllocator::size;

    //- Rewind the output buffer
    using OBufStreamAllocator::rewind;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
