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
    primaryMesh_(refCast<const fvMesh>(obr)),
    time_(primaryMesh_.time())
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
    
    // Set the actual output method here for the ADIOS group
    adios_select_method (groupID_, adiosMethod_.c_str(), methodParams_.c_str(), "");
    adios_allocate_buffer (ADIOS_BUFFER_ALLOC_NOW, 10);

    // Write initial conditions (including mesh)
    if (restartTime_ == 0.0) {
        write();
    } else {
        Info<< " Restart will happen later. We skip writing the initial 0 time " << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::adiosWrite::~adiosWrite()
{
    adios_free_group (groupID_); // not necessary but nice to cleanup
    adios_finalize (Pstream::myProcNo());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::adiosWrite::read_region(const dictionary& dict, regionInfo& r)
{
    Info<< " Region dict: " << r.name_ << endl;
    wordList toc_ (dict.toc());
    forAll(toc_, j)
    {
        Info<< "   TOC "<< j<<": " << toc_[j] << endl;
    }
    Info<< endl;
    dict.lookup("objectNames")  >> r.objectNames_;
    dict.lookup("cloudNames")   >> r.cloudNames_;
    dict.lookup("cloudAttribs") >> r.cloudAttribs_;

    // Set length of cell numbers array
    r.nCells_.setSize(Pstream::nProcs());
    // Set length of cell data size array
    r.cellDataSizes_.setSize(Pstream::nProcs());
    // Set length of particle numbers array
    r.nParticles_.setSize(Pstream::nProcs());

    // Do a basic check to see if the objectNames_ is accessible
    const fvMesh& m = time_.lookupObject<fvMesh>(r.name_);
    forAll(r.objectNames_, j)
    {
        if (m.foundObject<regIOobject>(r.objectNames_[j]))
        {
            Info<< " " << r.objectNames_[j];
        }
        else
        {
            WarningIn
                (
                 "Foam::adiosWrite::read(const dictionary&)"
                )   << "Object " << r.objectNames_[j] << " not found in "
                << "database. Available objects:" << nl << m.sortedToc()
                << endl;
        }

    }

    // Also print the cloud names
    forAll(r.cloudNames_, j)
    {
        Info<< " " << r.cloudNames_[j];
    }
    Info<< endl << endl;
}

void Foam::adiosWrite::read(const dictionary& dict)
{
    wordList toc_ (dict.toc());
    forAll(toc_, i)
    {
        Info<< " TOC "<< i<<": " << toc_[i] << endl;
    }
    Info<< endl << endl;

    /* Get the list of regions in adiosDict (as name list and dictionary list) */
    const entry* entryPtr = dict.lookupEntryPtr ("regions", false, false);
    if (entryPtr)
    {
        if (!entryPtr->isDict()) 
        {
            // a list of functionObjects
            PtrList<entry> r(entryPtr->stream());
            regions_.setSize(r.size());
            label nr = 0;
            forAllIter(PtrList<entry>, r, iter)
            {
                // safety:
                if (!iter().isDict())
                {
                    WarningIn ( "Foam::adiosWrite::read(const dictionary&)")   
                        << "Regions list should only contain dictionaries," 
                        << "one per region."
                        << endl;
                    continue;
                }
                //const word& key = iter().keyword();
                const dictionary& d = iter().dict();
                regions_[nr].name_ = iter().keyword();

                //regions_[nr].mesh_ = obr_.foundObject<fvMesh>(iter().keyword());
                //const fvMesh& m = time_.lookupObject<fvMesh>(iter().keyword());

                /* Process the dictionary here because we can't save the reference
                   into a list for later processing */
                /* Process each region, which should contain the fields and particles */

                read_region (d, regions_[nr]);
                nr++;
            }
        }
        else 
        {
            FatalIOErrorIn("Foam::adiosWrite::read(const dictionary&)", dict)
                << "List of regions should not be a dictionary but a list. "
                << exit(FatalIOError);
        }
    }
    else
    {
        FatalIOErrorIn("adiosWrite::read(const dictionary&)", dict)
            << "List of regions is required in adiosData in controlDict. "
            << exit(FatalIOError);
    }



    // Lookup in dictionary
    writeInterval_ = dict.lookupOrDefault<label>("writeInterval", 1);
    
    // Lookup chunk size if present
    adiosMethod_  = dict.lookupOrDefault<word>("adiosMethod", "MPI");
    methodParams_ = dict.lookupOrDefault<string>("methodparams", "");

    restartTime_ = dict.lookupOrDefault<scalar>("restartTime", 0.0);
    
    // Print info to terminal
    int writeprec = sizeof(ioScalar);
    Info<< type() << " " << name() << ":" << endl
        << "  Compiled with " << writeprec << " bytes precision." << endl
        << "  writing every " << writeInterval_ << " iterations:" << endl;
    if (restartTime_ > 0.0) 
    {
        Info<< "  restart time requested " << restartTime_ << endl;
    }
    
    
    // Check if writeInterval is a positive number
    if (writeInterval_ <= 0)
    {
        FatalIOErrorIn("adiosWrite::read(const dictionary&)", dict)
            << "Illegal value for writeInterval " << writeInterval_
            << ". It should be > 0."
            << exit(FatalIOError);
    }

    Info<< type() << " " << name() << ":" << endl
        << "  ADIOS output method: " << adiosMethod_ << endl 
        << "      with parameters: " << methodParams_ << endl
        << "       write interval: " << writeInterval_ << endl;
}


void Foam::adiosWrite::execute()
{
    /* execute() is called at the first non-zero time, after the calculation, before write() is called.
       This is a good point to do restart because fields created by other function objects exist at
       this point (e.g. fieldAverage variables)
    */
    static bool restarted = false; 
    Info<< "adiosWrite::execute() timeOutputValue  = " << obr_.time().timeOutputValue() << endl;
    if (restartTime_ > 0.0 && !restarted) 
    {
        Info<< "  restart time requested was " << restartTime_ << ". Let's do restart now." << endl;
        restarted = true; // even if it fails we don't try it again
        // classify fields in the object space, then read data from restart file for each
        classifyFields();
        if (!readData(restartTime_))
        {
            FatalErrorIn("adiosWrite::execute()")
                << "Restart reading failed for time " << restartTime_
                << exit(FatalIOError);
        }


        //runTime.setTime (restartTime_, 0 /* correct timeIndex is needed here */);
    }
}


void Foam::adiosWrite::end()
{
    // Nothing to be done here
    Info<< "adiosWrite::end() has been called at time " << obr_.time().timeName() << endl;
}


void Foam::adiosWrite::timeSet()
{
    // Nothing to be done here
    Info<< "adiosWrite::timeSet() has been called at time " << obr_.time().timeName() << endl;
}


void Foam::adiosWrite::defineVars()
{
    Info<< "adiosWrite::defineVars() has been called at time " << obr_.time().timeName() << endl;
    outputSize_ = 0;

    // define some info scalars: number of regions
    adios_define_var (groupID_, "nregions", "", adios_integer, NULL, NULL, NULL);
    outputSize_ += 4;

    forAll(regions_, regionID)
    {
        meshDefine(regionID); 
        fieldDefine(regionID);
        cloudDefine(regionID);
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
    Info<< "adiosWrite::write() has been called at time " << obr_.time().timeName() 
        << " time index " << obr_.time().timeIndex() << endl;
    // Check if we are going to write
    //if ( timeSteps_ == 0 )
    if ( timeSteps_ == nextWrite_ )
    {
        // Write info to terminal
        Info<< "Writing ADIOS data for time " << obr_.time().timeName() << endl;

        // Classify fields for all regions
        // at every write step in case new variables appear (they tend to, after the initial write at time 0)
        classifyFields();
    
#if 1
        if (timeSteps_ != 0)
        {
            // the size of patch variables is changing at every step, so in 
            // ADIOS we have to redefine the variables (all of them for simplicity)
            deleteDefinitions();
        }
        defineVars();
#else

        if (timeSteps_ == 0)
        {

            /* FIXME: code snippet to wait gdb attach on rank 0 process */
            /*
            Info<< "Pstream::myProcNo() = " << Pstream::myProcNo() << endl;
            if (Pstream::myProcNo() == 0) 
            {
                int z = 0;
                char hostname[256];
                gethostname(hostname, sizeof(hostname));
                int pid = getpid();
                WarningIn ("Foam::adiosWrite::write()")  << 
                    "PID " << pid << " on " << hostname << " ready for attach\n" << endl;
                while (0 == z)
                    sleep(5);
            }
            */
        
            // ADIOS requires to define all variables before writing anything
            Info<< "Define variables in ADIOS" << endl;
            defineVars();
        }
        else if (primaryMesh_.changing())
        {
            // Re-define all variables if mesh has changed (due to variable sizes change)
            Info<< "Redefine all variables in ADIOS because primary mesh has changed at time " << obr_.time().timeName() << endl;
            deleteDefinitions();
            defineVars();
        }
#endif

        // Create/reopen ADIOS output file, and tell ADIOS how man bytes we are going to write
        open();

        int n = regions_.size();
        adios_write (fileID_, "nregions", &n);

        forAll(regions_, regionID)
        {
            // Re-write mesh if dynamic or first time
            if (timeSteps_ == 0 || primaryMesh_.changing())
            {
                meshWrite(regionID); 
            }

            // Write field data
            fieldWrite(regionID);

            // Write cloud data 
            cloudWrite(regionID);
        }

        // Close the ADIOS dataset at every timestep. It must be closed!
        // Data is flushed from memory at this point
        close();
        
        // Calculate time of next write
        nextWrite_ = timeSteps_ + writeInterval_;
    }
    
    // Update counter
    timeSteps_++; 
}


// ************************************************************************* //
