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
    Foam::adiosCore

Description
    Core routines for reading/writing ADIOS files for OpenFOAM.

SourceFiles
    adiosCore.C
\*---------------------------------------------------------------------------*/

#ifndef adiosCore_H
#define adiosCore_H

#include "word.H"
#include "instantList.H"
#include "IOstreams.H"

#include "IBufStream.H"
#include "OBufStream.H"
#include "OCountStream.H"
#include "OCompactBufStream.H"
#include "OCompactCountStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class adiosCore Declaration
\*---------------------------------------------------------------------------*/

class adiosCore
{
private:

    // Private data

        //- Group name within ADIOS files
        word name_;

protected:

    // Protected data

        //- Directory name for ADIOS files
        static const word dataDirectory;

        //- File extension for ADIOS files
        static const word fileExt;

public:

    static const IOstream::streamFormat strFormat = IOstream::ASCII;


    // Constructors

        //- Construct for given ADIOS group name
        adiosCore(const word& groupName);


    //- Destructor
    virtual ~adiosCore();


    // Member Functions

        //- Return name of the ADIOS group
        const word& name() const
        {
            return name_;
        }


        //- Search a given directory for ADIOS time files
        static instantList findTimes
        (
            const fileName& directory = fileName::null,
            const word& constantName  = "constant"
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
