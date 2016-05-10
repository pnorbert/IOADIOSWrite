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
    Foam::OStreamAnalyze

Description
    An output stream with output to cerr

\*---------------------------------------------------------------------------*/

#ifndef OStreamAnalyze_H
#define OStreamAnalyze_H

#include "OSstream.H"
#include "OCompactStream.H"

#include <iostream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class reportstreambuf Declaration
\*---------------------------------------------------------------------------*/

//- A streambuf class for reporting bytes seen
class reportstreambuf
:
    public std::streambuf
{
    // Private data

    //- The number of bytes
    std::streamsize n_;

    const char* hexChars;


protected:
    // Protected members

        //- Put sequence of characters
        virtual std::streamsize xsputn(const char* s, std::streamsize num)
        {
            for (std::streamsize i=0; i < num; ++i)
            {
                std::cerr
                    << 'x' << (hexChars[((s[i] >> 4) & 0xF)])
                    << (hexChars[(s[i] & 0xF)]);
            }

            n_ += num;
            return num;
        }


        //- Set position pointer to absolute position
        //  For the counter, any positioning is ignored and it always acts like
        //  seekpos(0), which resets the count.
        virtual std::streampos seekpos
        (
            std::streampos sp,
            std::ios_base::openmode which = std::ios_base::in | std::ios_base::out
        )
        {
            n_ = 0;
            return 0;
        }


public:

    // Constructors

        //- Construct null
        reportstreambuf()
        :
            n_(0),
            hexChars("0123456789abcdef")
        {}


    // Access

        //- Get number of bytes counted
        std::streamsize size() const
        {
            return n_;
        }
};


/*---------------------------------------------------------------------------*\
                        Class reportstream Declaration
\*---------------------------------------------------------------------------*/

//- Trivial output stream for reporting bytes seen
class reportstream
:
    virtual public std::ios,
    public std::ostream
{
    // Private data

    reportstreambuf buf_;

public:

    // Constructors

        //- Construct null
        reportstream()
        :
            std::ostream(&buf_)
        {}


    // Member Functions

    // Access

        //- This hides both signatures of std::basic_ios::rdbuf()
        reportstreambuf* rdbuf()
        {
            return &buf_;
        }


        //- Get number of bytes counted
        std::streamsize size() const
        {
            return buf_.size();
        }


        //- Rewind the stream, reset the count
        void rewind()
        {
            buf_.pubseekpos(0);
            clear();    // for safety, clear any old errors
        }
};


/*---------------------------------------------------------------------------*\
                       Class OStreamAnalyze Declaration
\*---------------------------------------------------------------------------*/

//- An output stream for reporting bytes seen
class OStreamAnalyze
:
    public OSstream
{

    // Private Member Functions

        //- Disallow default bitwise copy construct
        OStreamAnalyze(const OCountStream&);

        //- Disallow default bitwise assignment
        void operator=(const OStreamAnalyze&);

        //- Base stream
        reportstream& baseStream()
        {
            return dynamic_cast<reportstream&>(stdStream());
        }


public:

    // Constructors

        //- Construct and set stream status
        OStreamAnalyze
        (
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            OSstream
            (
                *(new reportstream()),
                "OStreamAnalyze.sink",
                format,
                version
            )
        {}


    //- Destructor
    ~OStreamAnalyze()
    {
        delete &(baseStream());
    }


    // Member functions

    // Access

        //- Return the number of bytes counted
        std::streamsize size()
        {
            return baseStream().size();
        }


    // Edit

        //- Rewind the stream, reset the count, clearing any old errors
        void rewind()
        {
            baseStream().rewind();
            setGood();  // resynchronize with internal state
        }

};


/*---------------------------------------------------------------------------*\
                    Class OCompactStreamAnalyze Declaration
\*---------------------------------------------------------------------------*/

//- An output stream for reporting bytes seen, compact output
class OCompactStreamAnalyze
:
    public OCompactStream
{

    // Private Member Functions

        //- Disallow default bitwise copy construct
        OCompactStreamAnalyze(const OCompactStreamAnalyze&);

        //- Disallow default bitwise assignment
        void operator=(const OCompactStreamAnalyze&);

        //- Base stream
        reportstream& baseStream()
        {
            return dynamic_cast<reportstream&>(stdStream());
        }

public:

    // Constructors

        //- Construct and set stream status
        OCompactStreamAnalyze
        (
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            OCompactStream
            (
                *(new reportstream()),
                "OStreamAnalyze.sink",
                format,
                version
            )
        {}


    //- Destructor
    ~OCompactStreamAnalyze()
    {
        delete &(baseStream());
    }


    // Member functions

    // Access

        //- Return the number of bytes counted
        std::streamsize size()
        {
            return baseStream().size();
        }


    // Edit

        //- Rewind the stream, reset the count, clearing any old errors
        void rewind()
        {
            baseStream().rewind();
            setGood();  // resynchronize with internal state
        }

};


/*---------------------------------------------------------------------------*\
                      Class ORawStreamAnalyze Declaration
\*---------------------------------------------------------------------------*/

//- An output stream for reporting bytes seen, compact output
class ORawStreamAnalyze
:
    public OCompactStreamAnalyze
{
public:

    // Constructors

        //- Construct and set stream status
        ORawStreamAnalyze
        (
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            OCompactStreamAnalyze(format, version)
        {
            raw(true); // raw also forces binary format
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
