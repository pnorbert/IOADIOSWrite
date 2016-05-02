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
    Similar to OStringStream but using an externally managed buffer for
    its output.

    This allows the output buffer to be reused and can make it easier when
    writing out data.  It is the user's responsibility to ensure proper
    synchronization in the sizes. Provided that the external buffer is large
    enough that overflow does not occur, the following usage pattern
    should work.

    \code
        DynamicList<char> buffer(4096);     // allocate some large buffer

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
        DynamicList<char> buffer(4096);     // allocate some large buffer

        OBufStream os(buffer);
        os << "content1" << " and more content";
        buffer.setSize(os.size());          // synchronize sizes

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
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- A streambuf class for using externally allocated memory for its buffer
class memorybuf
:
    public std::streambuf
{
    // Private data

        //- Storage location
        char *buf_;

        //- The number of bytes in the storage
        std::streamsize len_;

protected:
    // Protected members

        //- Set position pointer to relative position
        virtual std::streampos seekoff
        (
            std::streamoff off,
            std::ios_base::seekdir way,
            std::ios_base::openmode which = std::ios_base::in | std::ios_base::out
        )
        {
            const bool testin  = which & std::ios_base::in;
            const bool testout = which & std::ios_base::out;

            if (way == std::ios_base::beg)
            {
                if (testin)
                {
                    setg(eback(), eback(), egptr());
                    gbump(off);
                }
                if (testout)
                {
                    setp(pbase(), epptr());
                    pbump(off);
                }

                return off;
            }

            if (way == std::ios_base::cur)
            {
                if (testin)
                {
                    gbump(off);
                }
                if (testout)
                {
                    pbump(off);
                }
            }
            else if (way == std::ios_base::end)
            {
                if (testin)
                {
                    gbump(off);
                }
                if (testout)
                {
                    pbump(off);
                }
            }

            if (testin)
            {
                return gptr() - eback();
            }
            if (testout)
            {
                return pptr() - pbase();
            }

            return -1;
        }


        //- Set position pointer to absolute position
        virtual std::streampos seekpos
        (
            std::streampos pos,
            std::ios_base::openmode which = std::ios_base::in | std::ios_base::out
        )
        {
            return seekoff(pos, std::ios_base::beg, which);
        }


        //- Get sequence of characters
        virtual std::streamsize xsgetn(char* s, std::streamsize n)
        {
            std::streamsize count = 0;

            // some optimization could be possible here
            while (count < n && gptr() < egptr())
            {
                *(s + count++) = *(gptr());
                gbump(1);
            }

            return count;
        }


        //- Put sequence of characters
        virtual std::streamsize xsputn(const char* s, std::streamsize n)
        {
            std::streamsize count = 0;

            // some optimization could be possible here
            while (count < n && pptr() < epptr())
            {
                *(pptr()) = *(s + count++);
                pbump(1);
            }

            return count;
        }


public:

    // Forward declarations of adapters

        class input;
        class output;


    // Constructors

        //- Construct for specified buffer
        memorybuf(char* buffer, std::streamsize num)
        :
            buf_(buffer),
            len_(num)
        {
            setg(buf_, buf_, buf_ + len_);
            setp(buf_, buf_ + len_);
        }

};


//- An allocation adapter for input
struct memorybuf::input : public memorybuf
{
    //- Construct with buffer and size
    input(const char *buffer, size_t nbytes)
    :
        memorybuf(const_cast<char*>(buffer), nbytes)
    {}


    //- Construct with buffer and size
    input(const void *buffer, size_t nbytes)
    :
        memorybuf
        (
            const_cast<char*>
            (
                reinterpret_cast<const char*>(buffer)
            ),
            nbytes
        )
    {}


    //- Construct with buffer and its inherent size
    template<class T>
    input(const UList<T>& buffer)
    :
        memorybuf
        (
            const_cast<char*>
            (
                reinterpret_cast<const char*>(buffer.cdata())
            ),
            buffer.size() * sizeof(T)
        )
    {}


    //- Construct with buffer and size
    template<class T>
    input(const UList<T>& buffer, label listSize)
    :
        memorybuf
        (
            const_cast<char*>
            (
                reinterpret_cast<const char*>(buffer.cdata())
            ),
            listSize * sizeof(T)
        )
    {}
};


//- An allocation adapter for output
struct memorybuf::output : public memorybuf
{
    //- Construct with buffer and size
    output(char *buffer, size_t nbytes)
    :
        memorybuf(buffer, nbytes)
    {}


    //- Construct with buffer and size
    output(void *buffer, size_t nbytes)
    :
        memorybuf
        (
            reinterpret_cast<char*>(buffer),
            nbytes
        )
    {}


    //- Construct with buffer and its inherent size
    template<class T>
    output(UList<T>& buffer)
    :
        memorybuf
        (
            reinterpret_cast<char*>(buffer.data()),
            buffer.size() * sizeof(T)
        )
    {}


    //- Construct with buffer and size
    template<class T>
    output(UList<T>& buffer, label listSize)
    :
        memorybuf
        (
            reinterpret_cast<char*>(buffer.data()),
            listSize * sizeof(T)
        )
    {}
};


/*---------------------------------------------------------------------------*\
                     Class OBufStreamAllocator Declaration
\*---------------------------------------------------------------------------*/

//- An allocator for std::ostream with external buffers
class OBufStreamAllocator
{
private:

    // Private data

        //- Reference to the underlying buffer
        memorybuf::output buf_;

protected:

    // Protected data

        //- The stream pointer
        std::ostream* strPtr_;

public:

    // Constructors

        //- Construct with buffer and number of bytes
        OBufStreamAllocator(char *buffer, size_t nbytes)
        :
            buf_(buffer, nbytes),
            strPtr_(new std::ostream(&buf_))
        {}


        //- Construct with buffer and number of bytes
        OBufStreamAllocator(void *buffer, size_t nbytes)
        :
            buf_(buffer, nbytes),
            strPtr_(new std::ostream(&buf_))
        {}


        //- Construct using specified buffer and its inherent storage size
        template<class T>
        OBufStreamAllocator(UList<T>& buffer)
        :
            buf_(buffer),
            strPtr_(new std::ostream(&buf_))
        {}


        //- Construct with buffer and size
        template<class T>
        OBufStreamAllocator(UList<T>& buffer, label listSize)
        :
            buf_(buffer, listSize),
            strPtr_(new std::ostream(&buf_))
        {}


    //- Destructor
    ~OBufStreamAllocator()
    {
        delete strPtr_;
    }


    // Member functions

        //- Return the current output position in the buffer
        std::streampos size()
        {
            return strPtr_->tellp();
        }


        //- Move to buffer start, clear errors
        void rewind()
        {
            strPtr_->rdbuf()->pubseekpos(0);
            strPtr_->clear(); // for safety, clear any old errors
        }

};


/*---------------------------------------------------------------------------*\
                         Class OBufStream Declaration
\*---------------------------------------------------------------------------*/

//- An OSstream with an attached external buffer
class OBufStream
:
    private OBufStreamAllocator,
    public OSstream
{
public:

    // Constructors

        //- Construct using specified buffer and number of bytes
        OBufStream
        (
            char* buffer,
            size_t nbytes,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            OBufStreamAllocator(buffer, nbytes),
            OSstream
            (
                *strPtr_,
                "OBufStream.sink",
                format,
                version
            )
        {}


        //- Construct using data area from buffer and number of bytes
        template<class T>
        OBufStream
        (
            UList<T>& buffer,
            size_t nbytes,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            OBufStreamAllocator(buffer.data(), nbytes),
            OSstream
            (
                *strPtr_,
                "OBufStream.sink",
                format,
                version
            )
        {}


        //- Construct using specified buffer and its inherent storage size
        template<class T>
        OBufStream
        (
            UList<T>& buffer,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            OBufStreamAllocator(buffer),
            OSstream
            (
                *strPtr_,
                "OBufStream.sink",
                format,
                version
            )
        {}


        //- Construct using specified buffer and its storage capacity
        template<class T>
        OBufStream
        (
            DynamicList<T>& buffer,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            OBufStreamAllocator(buffer, buffer.capacity()),
            OSstream
            (
                *strPtr_,
                "OBufStream.sink",
                format,
                version
            )
        {}


    //- Destructor
    ~OBufStream()
    {}


    // Member functions

        //- Return the current output position in the buffer
        using OBufStreamAllocator::size;

        //- Rewind the stream, clearing any old errors
        Ostream& rewind()
        {
            OBufStreamAllocator::rewind();
            setGood();  // resynchronize with internal state

            return *this;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
