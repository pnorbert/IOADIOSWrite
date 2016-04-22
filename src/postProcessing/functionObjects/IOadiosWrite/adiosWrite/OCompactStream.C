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

\*---------------------------------------------------------------------------*/

#include "OCompactStream.H"
#include "token.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::OCompactStream::OCompactStream
(
    ostream& os,
    const string& name,
    streamFormat format,
    versionNumber version,
    compressionType compression
)
:
    OSstream(os, name, format, version, compression),
    raw_(false)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


bool Foam::OCompactStream::raw() const
{
    return raw_ && format() != ASCII;
}


void Foam::OCompactStream::raw(bool b)
{
    raw_ = b;
    if (raw_)
    {
        format(IOstream::BINARY);
    }
}


Foam::Ostream& Foam::OCompactStream::write(const char c)
{
    if (format() == ASCII || c != token::NL)
    {
        // suppress newline for binary streams
        OSstream::write(c);
    }

    return *this;
}


Foam::Ostream& Foam::OCompactStream::write(const char* buf, std::streamsize count)
{
    if (format() != BINARY)
    {
        FatalIOErrorInFunction(*this)
            << "stream format not binary"
            << abort(FatalIOError);
    }

    if (raw_)
    {
        stdStream().write(buf, count);
        setState(stdStream().rdstate());
    }
    else
    {
        OSstream::write(buf, count);
    }

    return *this;
}


Foam::Ostream& Foam::OCompactStream::writeKeyword(const keyType& kw)
{
    write(kw);
    write(char(token::SPACE));

    return *this;
}


Foam::Ostream& Foam::OCompactStream::beginBlock(const word& keyword)
{
    write(keyword);
    beginBlock();

    return *this;
}


Foam::Ostream& Foam::OCompactStream::beginBlock()
{
    write(char(token::BEGIN_BLOCK));
    return *this;
}


Foam::Ostream& Foam::OCompactStream::endBlock()
{
    write(char(token::END_BLOCK));
    return *this;
}


void Foam::OCompactStream::indent()
{}


void Foam::OCompactStream::endl()
{
    write('\n');
    flush();
}


// ************************************************************************* //
