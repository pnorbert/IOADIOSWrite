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
    Foam::OCompactBufStream

Description
    A version of OStringStream that uses an externally managed buffer for
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

class OCompactBufStream
:
    private OBufStreamAllocator,
    public OCompactStream
{
public:

    // Constructors

        //- Construct using specified buffer for the output
        OCompactBufStream
        (
            const UList<char>& buffer,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            OBufStreamAllocator(buffer),
            OCompactStream
            (
                *strPtr_,
                "OStringStream.sinkFile", // keep same name as string-stream
                format,
                version
            )
        {}


        //- Construct using specified buffer for the output
        OCompactBufStream
        (
            const DynamicCharList& buffer,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            OBufStreamAllocator(buffer, buffer.capacity()),
            OCompactStream
            (
                *strPtr_,
                "OStringStream.sinkFile", // keep same name as string-stream
                format,
                version
            )
        {}


    //- Destructor
    ~OCompactBufStream()
    {}


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
