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
    Foam::OCompactCountStream

Description
    An output stream for calculating byte counts,
    using a more compact output (fewer spaces).
    Optionally supports raw binary output for data blocks.

\*---------------------------------------------------------------------------*/

#ifndef OCompactCountStream_H
#define OCompactCountStream_H

#include "OCompactStream.H"
#include "OCountStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class OCompactCountStream Declaration
\*---------------------------------------------------------------------------*/

class OCompactCountStream
:
    public OCompactStream
{

    // Private Member Functions

        //- Disallow default bitwise copy construct
        OCompactCountStream(const OCompactCountStream&);

        //- Disallow default bitwise assignment
        void operator=(const OCompactCountStream&);

public:

    // Constructors

        //- Construct and set stream status
        OCompactCountStream
        (
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            OCompactStream
            (
                *(new ocountstream()),
                "OCompactCountStream.sink",
                format,
                version
            )
        {}


    //- Destructor
    ~OCompactCountStream()
    {
        delete &dynamic_cast<ocountstream&>(stdStream());
    }


    // Member functions

    // Access

        //- Return the number of bytes counted
        std::streamsize size()
        {
            return dynamic_cast<ocountstream&>(stdStream()).size();
        }


    // Edit

        //- Rewind the stream, reset the count, clearing any old errors
        void rewind()
        {
            dynamic_cast<ocountstream&>(stdStream()).rewind();
            setGood();  // resynchronize with internal state
        }

};


/*---------------------------------------------------------------------------*\
                       Class ORawCountStream Declaration
\*---------------------------------------------------------------------------*/

//- An output stream for calculating byte counts,
//  using a more compact output (fewer spaces).
//  Supports raw binary output for data blocks.
class ORawCountStream
:
    public OCompactCountStream
{
public:

    // Constructors

        //- Construct and set stream status
        ORawCountStream
        (
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            OCompactCountStream(format, version)
        {
            raw(true); // raw also forces binary format
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
