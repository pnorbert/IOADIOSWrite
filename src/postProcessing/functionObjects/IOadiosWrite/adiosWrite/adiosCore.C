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

\*---------------------------------------------------------------------------*/

#include "adiosCore.H"

#include "fileNameList.H"
#include "IStringStream.H"
#include "token.H"
#include "OSspecific.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word
Foam::adiosCore::dataDirectory("adiosData");

const Foam::word
Foam::adiosCore::fileExt("bp");

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adiosCore::adiosCore
(
    const word& name
)
:
    name_(name)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::adiosCore::~adiosCore()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::instantList Foam::adiosCore::findTimes
(
    const fileName& directory,
    const word& constantName
)
{
    // Read directory entries into a list
    fileNameList dirEntries(readDir(directory/dataDirectory, fileName::FILE, false));

    // Initialise instant list
    instantList Times(dirEntries.size() + 1);
    label nTimes = 0;

    // Check for "constant" - not yet useful
    bool haveConstant = false;
    forAll(dirEntries, i)
    {
        if
        (
            dirEntries[i].ext() == fileExt
         && dirEntries[i].lessExt() == constantName
        )
        {
            Times[nTimes].value() = 0;
            Times[nTimes].name()  = dataDirectory/dirEntries[i];
            nTimes++;
            haveConstant = true;
            break;
        }
    }

    // Read and parse all the entries in the directory
    forAll(dirEntries, i)
    {
        if
        (
            dirEntries[i].ext() == fileExt
        )
        {
            IStringStream timeStream(dirEntries[i].lessExt());
            token timeToken(timeStream);

            if (timeToken.isNumber() && timeStream.eof())
            {
                Times[nTimes].value() = timeToken.number();
                Times[nTimes].name()  = dataDirectory/dirEntries[i];
                nTimes++;
            }
        }
    }

    // Reset the length of the times list
    Times.setSize(nTimes);

    if (haveConstant)
    {
        if (nTimes > 2)
        {
            std::sort(&Times[1], Times.end(), instant::less());
        }
    }
    else if (nTimes > 1)
    {
        std::sort(&Times[0], Times.end(), instant::less());
    }

    return Times;
}


/// Foam::word Foam::Time::findInstancePath(const instant& t) const
/// {
///     const fileName directory = path();
///     const word& constantName = constant();
///
///     // Read directory entries into a list
///     fileNameList dirEntries(readDir(directory, fileName::DIRECTORY));
///
///     forAll(dirEntries, i)
///     {
///         scalar timeValue;
///         if (readScalar(dirEntries[i].c_str(), timeValue) && t.equal(timeValue))
///         {
///             return dirEntries[i];
///         }
///     }
///
///     if (t.equal(0.0))
///     {
///         // Looking for 0 or constant. 0 already checked above.
///         if (isDir(directory/constantName))
///         {
///             return constantName;
///         }
///     }
///
///     return word::null;
/// }

// ************************************************************************* //
