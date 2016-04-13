/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |               2015 Norbert Podhorszki
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

#include "adiosWrite.H"
#include "fileName.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::adiosWrite::open()
{
    Info<< "adiosWrite::open:" << endl;

    // Create output directory if nonexisting in the beginning
    static bool checkdir = true;
    if (checkdir && !isDir(dataDirectory))
    {
        mkDir(dataDirectory);
    }
    checkdir = false;


    // Find and create filename
    char mode[] = "w";

    fileName dataFile(dataDirectory/obr_.time().timeName() + ".bp");

    // Create/open a new file collectively.
    /*
    static int i = 0;
    if (timeSteps_ == 0)
    {
        mode[0] = 'w';
        do
        {
            //sprintf(dataFile, "adiosData/%s%i.bp", name_.c_str(), i);
            sprintf(dataFile, "%s%i.bp", name_.c_str(), i);
            i++;
        }
        while (isFile(dataFile));
    }
    else
    {
        mode[0] = 'a';
        sprintf(dataFile, "%s%i.bp", name_.c_str(), i);
    }
    */

    // Print info to terminal
    Info<< "  adiosWrite: Chosen filename " << dataFile << endl << endl;

    int err = adios_open
    (
        &fileID_,
        name().c_str(),     // group-name
        dataFile.c_str(),
        mode,
        comm_
    );

    if (!err)
    {
        // Tell ADIOS how many bytes we are going to write in this step (by this process)
        uint64_t total_size; // user data bytes + metadata, buffer size should be bigger
        adios_group_size(fileID_, outputSize_, &total_size);
    }
    else
    {
        WarningInFunction
            << "File " << name()
            << " could not be created: " << nl
            << adios_get_last_errmsg();
        fileID_ = 0;
    }
}


void Foam::adiosWrite::close()
{
    MPI_Barrier(MPI_COMM_WORLD);
    Info<< "adiosWrite::close" << endl << endl;

    // Close the file
    if (fileID_)
    {
        adios_close(fileID_);
    }
}



// ************************************************************************* //
