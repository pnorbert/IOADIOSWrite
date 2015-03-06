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
#include "dictionary.H"
#include "scalar.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(adiosWrite, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adiosWrite::adiosWrite
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    mesh_(refCast<const fvMesh>(obr))
{
    Info<< "adiosWrite constructor called " << endl;
    MPI_Comm_dup (MPI_COMM_WORLD, &comm_);

    // Initilize ADIOS
    adios_init_noxml (comm_);
    // Create a group to hold all variable definitions
    adios_declare_group (&groupID_, name_.c_str(), "", adios_flag_yes);

    // Set next write NOW
    nextWrite_ = 0;
    timeSteps_ = 0;
    
    // Read dictionary
    read(dict);
    
    // Classify fields
    nFields_ = classifyFields();
    
    // Set length of cell numbers array
    nCells_.setSize(Pstream::nProcs());
    // Set length of cell data size array
    cellDataSizes_.setSize(Pstream::nProcs());

    // Set the actual output method here for the ADIOS group
    adios_select_method (groupID_, adiosMethod_.c_str(), methodParams_.c_str(), "");
    adios_allocate_buffer (ADIOS_BUFFER_ALLOC_NOW, 10);

    // Only do if some fields are to be written
    if (nFields_ || cloudNames_.size())
    {
        // Write initial conditions (including mesh)
        write();
    }


}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::adiosWrite::~adiosWrite()
{
    adios_free_group (groupID_); // not necessary but nice to cleanup
    adios_finalize (Pstream::myProcNo());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::adiosWrite::read(const dictionary& dict)
{
    // Lookup in dictionary
    dict.lookup("objectNames") >> objectNames_;
    dict.lookup("cloudNames") >> cloudNames_;
    dict.lookup("cloudAttribs") >> cloudAttribs_;
    dict.lookup("writeInterval") >> writeInterval_;
    
    // Lookup chunk size if present
    adiosMethod_  = dict.lookupOrDefault<word>("adiosMethod", "MPI");
    methodParams_ = dict.lookupOrDefault<string>("methodparams", "");
    
    // Print info to terminal
    int writeprec = sizeof(ioScalar);
    Info<< type() << " " << name() << ":" << endl
        << "  Compiled with " << writeprec << " bytes precision." << endl
        << "  writing every " << writeInterval_ << " iterations:"
        << endl
        << "   ";
    
    // Do a basic check to see if the objectNames_ is accessible
    forAll(objectNames_, i)
    {
        if (obr_.foundObject<regIOobject>(objectNames_[i]))
        {
            Info<< " " << objectNames_[i];
        }
        else
        {
            WarningIn
            (
                "Foam::writeRegisteredObject::read(const dictionary&)"
            )   << "Object " << objectNames_[i] << " not found in "
                << "database. Available objects:" << nl << obr_.sortedToc()
                << endl;
        }

    }
    
    // Also print the cloud names
    forAll(cloudNames_, i)
    {
      Info<< " " << cloudNames_[i];
    }
    Info<< endl << endl;
    
    // Check if writeInterval is a positive number
    if (writeInterval_ <= 0)
    {
        FatalIOErrorIn("adiosWrite::read(const dictionary&)", dict)
            << "Illegal value for writeInterval " << writeInterval_
            << ". It should be > 0."
            << exit(FatalIOError);
    }
}


void Foam::adiosWrite::execute()
{
    // Nothing to be done here
    Info<< "adiosWrite::execute() has been called " << endl;
}


void Foam::adiosWrite::end()
{
    // Nothing to be done here
    Info<< "adiosWrite::end() has been called " << endl;
}


void Foam::adiosWrite::timeSet()
{
    // Nothing to be done here
    Info<< "adiosWrite::timeSet() has been called " << endl;
}


void Foam::adiosWrite::defineVars()
{
    Info<< "adiosWrite::defineVars() has been called at time " << obr_.time().timeName() << endl;
    outputSize_ = 0;
    if (nFields_)
    {
        meshDefine();
        fieldDefine();
    }
    if (cloudNames_.size()) {
        cloudDefine();
    }
}

void Foam::adiosWrite::deleteDefinitions()
{
    Info<< "adiosWrite::deleteDefinitions() has been called at time " << obr_.time().timeName() << endl;
    // In ADIOS we need to remove all variable definitions in order
    // to make a new list of definitions in case the mesh changes
    adios_delete_vardefs (groupID_);
}

void Foam::adiosWrite::write()
{
    Info<< "adiosWrite::write() has been called at time " << obr_.time().timeName() << endl;
    // Check if we are going to write
    if ( timeSteps_ == 0 )
    //if ( timeSteps_ == nextWrite_ )
    {
        // Write info to terminal
        Info<< "Writing ADIOS data for time " << obr_.time().timeName() << endl;

        if (timeSteps_ == 0)
        {
            // ADIOS requires to define all variables before writing anything
            Info<< "Define variables in ADIOS" << endl;
            defineVars();
        }
        else if (mesh_.changing())
        {
            // Re-define all variables if mesh has changed (due to variable sizes change)
            Info<< "Redefine all variables in ADIOS because mesh has changed at time " << obr_.time().timeName() << endl;
            deleteDefinitions();
            defineVars();
        }

        // Create/reopen ADIOS output file
        open();

        // Only write field data if fields are specified
        if (nFields_)
        {
            // Re-write mesh if dynamic or first time
            if (timeSteps_ == 0 || mesh_.changing())
            {
                meshWrite();
            }
            
            // Write field data
            //fieldWrite();
        }
        
        // Only write cloud data if any clouds is specified
        if (cloudNames_.size())
        {
            cloudWrite();
        }
        
        // Close the ADIOS dataset at every timestep
        // Data is flushed from memory at this point
        close();
        
        // Calculate time of next write
        nextWrite_ = timeSteps_ + writeInterval_;
    }
    
    // Update counter
    timeSteps_++; 
}


// ************************************************************************* //
