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
    Foam::OCompactStream

Description
    A generic output stream using an output with reduced spaces

SourceFiles
    OCompactStream.C

\*---------------------------------------------------------------------------*/

#ifndef OCompactStream_H
#define OCompactStream_H

#include "OSstream.H"
#include "fileName.H"
#include <iostream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class OCompactStream Declaration
\*---------------------------------------------------------------------------*/

class OCompactStream
:
    public OSstream
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        OCompactStream(const OCompactStream&);

        //- Disallow default bitwise assignment
        void operator=(const OCompactStream&);

public:

    // Constructors

        //- Set stream status
        OCompactStream
        (
            ostream& os,
            const string& name,
            streamFormat format=ASCII,
            versionNumber version=currentVersion,
            compressionType compression=UNCOMPRESSED
        );


    // Member functions

        // Write functions

            //- Add indentation characters
            virtual void indent();

            //- Write the keyword followed by an appropriate indentation
            virtual Ostream& writeKeyword(const keyType&);

            //- Write begin block group with the given name
            //  Uses the appropriate indentation,
            //  does not include a trailing newline.
            virtual Ostream& beginBlock(const word&);

            //- Write begin block group without a name
            //  Uses the appropriate indentation,
            //  does not include a trailing newline.
            virtual Ostream& beginBlock();

            //- Write end block group
            //  Uses the appropriate indentation,
            //  does not include a trailing newline.
            virtual Ostream& endBlock();


        // Stream state functions

            //- Add newline and flush stream
            virtual void endl();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
