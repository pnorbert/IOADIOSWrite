/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd
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
    Foam::IStringStreamBuf

Description
    A non-copying version of IStringStream.
    This allows the input buffer to be filled and refilled from
    various sources

Notes
    The buffer is implemented as a DynamicList instead of a C++ string to
    ensure that the data are continguous in memory.
    In C++11, a string is stored continguously in memory, but includes a
    terminating null-character as well.
\*---------------------------------------------------------------------------*/

#ifndef IStringStreamBuf_H
#define IStringStreamBuf_H

#include "IStringStream.H"
#include "DynamicCharList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class IStringStreamBuf Declaration
\*---------------------------------------------------------------------------*/

class IStringStreamBuf
:
    public IStringStream,
    private DynamicCharList
{
private:

    typedef DynamicCharList StorageList;

    //- Attach storage buffer to the stream and allow nbytes of data
    //  When nbytes=0, allow usage of the entire storage buffer
    //  Rewinds the buffer for internal consistency
    void resized(size_t nbytes)
    {
        if (nbytes)
        {
            stdStream().rdbuf()->pubsetbuf(data(), nbytes);
        }
        else
        {
            stdStream().rdbuf()->pubsetbuf(data(), capacity());
        }

        stdStream().rdbuf()->pubseekpos(0); // rewind
        stdStream().clear();    // clear old errors
        setGood();              // resynchronize with internal state
        lineNumber(0);
    }

public:

    // Constructors

        //- Construct null, use 'reserve' to specify the desired buffer size
        IStringStreamBuf
        (
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            IStringStream("", format, version), // an empty buffer to avoid copy-in behaviour
            StorageList()
        {
            reserve();
        }


        //- Construct with the specified initial buffer size
        IStringStreamBuf
        (
            size_t nbytes,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            IStringStream("", format, version), // an empty buffer to avoid copy-in behaviour
            StorageList()
        {
            reserve(nbytes);
        }


    //- Destructor
    ~IStringStreamBuf()
    {}


    // Member functions

        // Access

            //- Read access to the underlying storage
            using StorageList::cdata;

            //- Read/write access to the underlying storage
            using StorageList::data;

            //- Size of the underlying storage.
            using StorageList::capacity;


            //- Clear the list and delete storage.
            void clearStorage()
            {
                StorageList::clearStorage();
                reserve();
            }


            //- Reserve allocation space for at least this size
            //  Never shrinks the allocated size.
            //  Attach storage buffer to the stream, automatically rewinds the buffer
            void reserve(const label nbytes = 0)
            {
                StorageList::reserve(nbytes);
                resized(nbytes);
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
