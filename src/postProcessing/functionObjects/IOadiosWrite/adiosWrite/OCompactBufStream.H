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
    Foam::OCompactBufStream

Description
    Similar to OStringStream but using an externally managed buffer for
    its output and using a more compact output (fewer spaces)

\*---------------------------------------------------------------------------*/

#ifndef OCompactBufStream_H
#define OCompactBufStream_H

#include "OBufStream.H"
#include "OCompactStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class OCompactBufStream Declaration
\*---------------------------------------------------------------------------*/

//- An OSstream with an attached external buffer and compact output
class OCompactBufStream
:
    private OBufStreamAllocator,
    public OCompactStream
{
public:

    // Constructors

        //- Construct using specified buffer and number of bytes
        OCompactBufStream
        (
            char* buffer,
            size_t nbytes,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            OBufStreamAllocator(buffer, nbytes),
            OCompactStream
            (
                *strPtr_,
                "OBufStream.sink",
                format,
                version
            )
        {}


        //- Construct using data area from buffer and number of bytes
        template<class T>
        OCompactBufStream
        (
            UList<T>& buffer,
            size_t nbytes,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            OBufStreamAllocator(buffer.data(), nbytes),
            OCompactStream
            (
                *strPtr_,
                "OBufStream.sink",
                format,
                version
            )
        {}


        //- Construct using specified buffer and its inherent storage size
        template<class T>
        OCompactBufStream
        (
            UList<T>& buffer,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            OBufStreamAllocator(buffer),
            OCompactStream
            (
                *strPtr_,
                "OBufStream.sink",
                format,
                version
            )
        {}


        //- Construct using specified buffer and its storage capacity
        template<class T>
        OCompactBufStream
        (
            DynamicList<T>& buffer,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            OBufStreamAllocator(buffer, buffer.capacity()),
            OCompactStream
            (
                *strPtr_,
                "OBufStream.sink",
                format,
                version
            )
        {}


    //- Destructor
    ~OCompactBufStream()
    {}


    // Member functions

        //- Return the current output position in the buffer
        using OBufStreamAllocator::size;

        //- Rewind the output buffer
        using OBufStreamAllocator::rewind;

};


/*---------------------------------------------------------------------------*\
                        Class ORawBufStream Declaration
\*---------------------------------------------------------------------------*/

//- An OSstream with an attached external buffer using a more compact output
//  (fewer spaces) and raw binary output
class ORawBufStream
:
    public OCompactBufStream
{
public:

    // Constructors

        //- Construct using specified buffer and number of bytes
        ORawBufStream
        (
            char* buffer,
            size_t nbytes,
            streamFormat format=BINARY,
            versionNumber version=currentVersion
        )
        :
            OCompactBufStream(buffer, nbytes, format, version)
        {
            raw(true);
        }


        //- Construct using data area from buffer and number of bytes
        template<class T>
        ORawBufStream
        (
            UList<T>& buffer,
            size_t nbytes,
            streamFormat format=BINARY,
            versionNumber version=currentVersion
        )
        :
            OCompactBufStream(buffer, nbytes, format, version)
        {
            raw(true);
        }


        //- Construct using specified buffer and its inherent storage size
        template<class T>
        ORawBufStream
        (
            UList<T>& buffer,
            streamFormat format=BINARY,
            versionNumber version=currentVersion
        )
        :
            OCompactBufStream(buffer, format, version)
        {
            raw(true);
        }


        //- Construct using specified buffer and its storage capacity
        template<class T>
        ORawBufStream
        (
            DynamicList<T>& buffer,
            streamFormat format=BINARY,
            versionNumber version=currentVersion
        )
        :
            OCompactBufStream(buffer, format, version)
        {
            raw(true);
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
