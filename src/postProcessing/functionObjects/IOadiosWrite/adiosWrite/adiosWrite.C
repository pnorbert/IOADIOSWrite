/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |               2015 Norbert Podhorszki
                            |               2016 OpenCFD Ltd.
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
#include "demandDrivenData.H"

// some internal pre-processor stringifications
#undef STRINGIFY
#undef TO_STRING
#define STRINGIFY(x) #x
#define TO_STRING(x) STRINGIFY(x)

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(adiosWrite, 0);
}


void Foam::adiosWrite::read_region(const dictionary& dict, regionInfo& rInfo)
{
    Info<< " Region dict: " << rInfo.name_ << endl;
    wordList toc(dict.toc());
    forAll(toc, j)
    {
        Info<< "   TOC " << j << ": " << toc[j] << endl;
    }
    Info<< endl;
    dict.lookup("objectNames")  >> rInfo.objectNames_;
    dict.lookup("cloudNames")   >> rInfo.cloudNames_;
    dict.lookup("cloudAttribs") >> rInfo.cloudAttribs_;

#ifdef FOAM_ADIOS_CELL_SHAPES
    // Set length of cell data size array
    rInfo.cellDataSizes_.setSize(Pstream::nProcs());
#endif

    // Set length of particle numbers array
    rInfo.nParticles_.setSize(Pstream::nProcs());

    // Do a basic check to see if the objectNames_ is accessible
    const fvMesh& mesh = time_.lookupObject<fvMesh>(rInfo.name_);

    DynamicList<word> missingObjects(rInfo.objectNames_.size());
    forAll(rInfo.objectNames_, i)
    {
        if (mesh.foundObject<regIOobject>(rInfo.objectNames_[i]))
        {
            Info<< " " << rInfo.objectNames_[i];
        }
        else
        {
            missingObjects.append(rInfo.objectNames_[i]);
        }
    }

    if (missingObjects.size())
    {
        WarningInFunction
            << "Objects not found in database:" << missingObjects
            << "Available objects:" << nl << mesh.sortedToc()
            << endl;
    }

    // Also print the cloud names
    forAll(rInfo.cloudNames_, i)
    {
        Info<< " " << rInfo.cloudNames_[i];
    }
    Info<< nl << endl;
}


size_t Foam::adiosWrite::defineVars(bool updateMesh)
{
    Info<< "adiosWrite::defineVars(" << (updateMesh ? "updateMesh" : "")
        << ") has been called at time "
        << obr_.time().timeName()
        << " time index " << obr_.time().timeIndex() << endl;

    outputSize_ = 0;
    size_t maxLen = 0;
    size_t bufLen = 0;

    // OpenFOAM build information (also contains version)
    defineAttribute
    (
        "version",
        adiosCore::foamAttribute,
        Foam::FOAMbuild
    );

    // OpenFOAM platform tag (WM_ARCH + WM_COMPILER)
    defineAttribute
    (
        "platform",
        adiosCore::foamAttribute,
#ifdef FOAM_PLATFORM
        TO_STRING(FOAM_PLATFORM)
#else
        ""
# error "FOAM_PLATFORM not defined"
#endif
    );

    // OpenFOAM label size (32|64)
    defineIntAttribute
    (
        "label-size",
        adiosCore::foamAttribute,
        adiosTraits<label>::nBits
    );

    // OpenFOAM scalar type (single-precision?)
    defineAttribute
    (
        "precision",
        adiosCore::foamAttribute,
        adiosTraits<scalar>::precisionName
    );

    // mesh updated - may also want moved points etc.
//  enum readUpdateState
//  {
//     UNCHANGED,
//     POINTS_MOVED,
//     TOPO_CHANGE,
//     TOPO_PATCH_CHANGE
//  };

    // store as integer instead of bool: minimal overhead, more flexibility
    defineIntAttribute
    (
        "updateMesh",
        adiosCore::foamAttribute,
        updateMesh
    );

    // other general information
    defineIntAttribute("nProcs",    adiosCore::foamAttribute, Pstream::nProcs());
    defineIntAttribute("nregions",  adiosCore::foamAttribute, regions_.size());

    // General information (as variable)
    defineIntVariable("nregions");      // number of regions
    defineIntVariable("time/index");
    defineScalarVariable("time/value");

    forAll(regions_, regionI)
    {
        regionInfo& rInfo = regions_[regionI];
        const fileName varPath = "region" + Foam::name(rInfo.index_);

        const fvMesh& mesh = time_.lookupObject<fvMesh>(rInfo.name_);
        const polyPatchList& patches = mesh.boundaryMesh();

        defineAttribute("name", varPath, rInfo.name_);
        defineIntAttribute("npatches", varPath, patches.size());

        forAll(patches, patchI)
        {
            const polyPatch& p = patches[patchI];
            fileName patchPath = varPath/"patch" + Foam::name(patchI);

            // global attributes for this patch
            defineAttribute("name",  patchPath, p.name());
            defineAttribute("type",  patchPath, p.type());
        }

        if (updateMesh)
        {
            bufLen = meshDefine(rInfo);
            maxLen = Foam::max(maxLen, bufLen);
        }

        bufLen = fieldDefine(rInfo);
        maxLen = Foam::max(maxLen, bufLen);

        bufLen = cloudDefine(rInfo);
        maxLen = Foam::max(maxLen, bufLen);
    }

    return maxLen;
}


void Foam::adiosWrite::deleteDefinitions()
{
    Info<< "adiosWrite::deleteDefinitions() has been called at time "
        << obr_.time().timeName()
        << " time index " << obr_.time().timeIndex() << endl;

    // In ADIOS we need to remove all variable definitions in order
    // to make a new list of definitions in case the mesh changes
    adios_delete_vardefs(groupID_);

    // also cleanup all old attributes
    adios_delete_attrdefs(groupID_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adiosWrite::adiosWrite
(
    const word& groupName,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    adiosCore(groupName),
    shapeLookupPtr_(NULL),
    obr_(obr),
    primaryMesh_(refCast<const fvMesh>(obr)),
    time_(primaryMesh_.time())
{
    Info<< "adiosWrite constructor called " << endl;
    MPI_Comm_dup(MPI_COMM_WORLD, &comm_);

    // Initilize ADIOS
    adios_init_noxml(comm_);

    // Create a group to hold all variable definitions
    adios_declare_group(&groupID_, name().c_str(), "", adios_flag_yes);

    // Set next write NOW
    nextWrite_ = 0;
    timeSteps_ = 0;

    // Read dictionary
    read(dict);

    // Set the actual output method here for the ADIOS group
    adios_select_method
    (
        groupID_,
        adiosMethod_.c_str(),
        methodParams_.c_str(),
        ""
    );
    adios_allocate_buffer(ADIOS_BUFFER_ALLOC_NOW, 10);

    // Write initial conditions (including mesh) if restartTime not set
    if (restartTime_ == VGREAT)
    {
        write();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::adiosWrite::~adiosWrite()
{
    adios_free_group(groupID_); // not necessary but nice to cleanup
    adios_finalize(Pstream::myProcNo());

    deleteDemandDrivenData(shapeLookupPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::adiosWrite::read(const dictionary& dict)
{
    wordList toc_ (dict.toc());
    forAll(toc_, i)
    {
        Info<< " TOC "<< i<<": " << toc_[i] << endl;
    }
    Info<< endl << endl;

    // Get the list of regions in adiosDict (as name list and dictionary list)
    const entry* entryPtr = dict.lookupEntryPtr("regions", false, false);

    if (entryPtr->isDict())
    {
        // Regions dictionary
        const dictionary& allRegionsDict = entryPtr->dict();

        label nRegion = 0;
        regions_.setSize(allRegionsDict.toc().size());
        forAllConstIter(dictionary, allRegionsDict, iter)
        {
            if (!iter().isDict())
            {
                FatalIOErrorInFunction(allRegionsDict)
                    << "Regions must be specified in dictionary format with "
                    << "the region name as the keyword"
                    << exit(FatalIOError);
            }

            const dictionary& regionDict = iter().dict();
            regions_[nRegion].index_ = nRegion;
            regions_[nRegion].name_  = iter().keyword();

            // Process each region, which should contain the fields and particles
            read_region(regionDict, regions_[nRegion]);
            ++nRegion;
        }
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "regions must be specified in dictionary format"
            << exit(FatalIOError);
    }

    // Lookup in dictionary
    writeInterval_ = dict.lookupOrDefault<label>("writeInterval", 1);

    // Lookup chunk size if present
    adiosMethod_  = dict.lookupOrDefault<word>("adiosMethod", "MPI");
    methodParams_ = dict.lookupOrDefault<string>("methodparams", "");

    // Print info to terminal
    Info<< type() << " " << name() << ":" << endl
        << "  Compiled with " << adiosTraits<scalar>::adiosSize << " bytes precision." << endl;

    if (sizeof(ioScalar) != adiosTraits<scalar>::adiosSize)
    {
        Info<< "  Writing IO with " << sizeof(ioScalar) << " bytes precision." << endl;
    }
    Info<< "  writing every " << writeInterval_ << " iterations:" << endl;


    restartTime_ = VGREAT;
    if (dict.readIfPresent("restartTime", restartTime_))
    {
        Info<< "  Restart time requested " << restartTime_ << endl;
    }


    // Check if writeInterval is a positive number
    if (writeInterval_ <= 0)
    {
        FatalIOErrorInFunction(dict)
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
    // execute() is called at the first non-zero time, after the calculation,
    // before write() is called.
    // This is a good point to do restart because fields created by other
    // function objects exist at this point (e.g. fieldAverage variables)

    static bool restarted = false;

    Info<< "adiosWrite::execute() timeOutputValue  = "
        << obr_.time().timeOutputValue() << endl;

    const scalar dt = obr_.time().deltaTValue();
    if (!restarted && (obr_.time().value() + 0.5*dt > restartTime_))
    {
        Info<< "  restart time requested was " << restartTime_
            << ". Let's do restart now." << endl;

        restarted = true; // even if it fails we don't try it again

        // Classify fields in the object space, then read data from restart
        // file for each
        classifyFields();

        Time& time = const_cast<Time&>(obr_.time());
        const instantList times = time.times();
        const label timeIndex = Time::findClosestTimeIndex(times, restartTime_);
        time.setTime(times[timeIndex], timeIndex);

        Info<<"times: " << times << endl;

        Info<<"look for times in " << dataDirectory << endl;
        const instantList adiosTimes = adiosCore::findTimes();

        Info<<"adios-times: " << adiosTimes << endl;
        const label adiosTimeIndex = Time::findClosestTimeIndex(adiosTimes, restartTime_);
        if (adiosTimeIndex < 0)
        {
            FatalErrorInFunction
                << "No appropriate ADIOS restart found for time " << restartTime_
                << exit(FatalIOError);
        }
        else
        {
            Info<<"restart from adios " << adiosTimes[adiosTimeIndex] << endl;
            if (!readData(adiosTimes[adiosTimeIndex]))
            {
                FatalErrorInFunction
                    << "Restart reading failed for time " << restartTime_
                    << exit(FatalIOError);
            }
        }
    }
}


void Foam::adiosWrite::end()
{
    // Nothing to be done here
    Info<< "adiosWrite::end() has been called at time "
        << obr_.time().timeName()
        << " time index " << obr_.time().timeIndex() << endl;
}


void Foam::adiosWrite::timeSet()
{
    // Nothing to be done here
    Info<< "adiosWrite::timeSet() has been called at time "
        << obr_.time().timeName()
        << " time index " << obr_.time().timeIndex() << endl;
}


void Foam::adiosWrite::write()
{
    Info<< "adiosWrite::write() has been called at "
        << "time " << obr_.time().timeName()
        << " time index " << obr_.time().timeIndex() << endl;

    size_t maxLen = 0;
    size_t bufLen = 0;

    // Check if we are going to write
    //if ( timeSteps_ == 0 )
    if (timeSteps_ == nextWrite_)
    {
        // Write info to terminal
        Info<< "Writing ADIOS data for time " << obr_.time().timeName() << endl;

        // Re-write mesh if dynamic or first time
        const bool updateMesh = (timeSteps_ == 0 || primaryMesh_.changing());

        // Classify fields for all regions at every write step in case new
        // variables appear, e.g. via function objects
        classifyFields();

#if 1
        if (timeSteps_ != 0)
        {
            // The size of patch variables is changing at every step, so in
            // ADIOS we have to redefine the variables
            // Note: updating all of them for simplicity
            deleteDefinitions();
        }

        bufLen = defineVars(updateMesh);
        maxLen = Foam::max(maxLen, bufLen);
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
            bufLen = defineVars(updateMesh);
            maxLen = Foam::max(maxLen, bufLen);
        }
        else if (primaryMesh_.changing())
        {
            // Re-define all variables if mesh has changed
            Info<< "Redefine all variables in ADIOS because primary mesh has "
                << "changed at time " << obr_.time().timeName() << endl;
            deleteDefinitions();
            bufLen = defineVars(updateMesh);
            maxLen = Foam::max(maxLen, bufLen);
        }
#endif

        // Create/reopen ADIOS output file, and tell ADIOS how many bytes we
        // are going to write
        open();

        Pout<<"reserve write buffer " << maxLen << endl;
        iobuffer_.reserve(maxLen);

        // General information (as variable)
        writeIntVariable("nregions", regions_.size());
        writeIntVariable("time/index",    time_.timeIndex());
        writeScalarVariable("time/value", time_.timeOutputValue());

        forAll(regions_, regionI)
        {
            regionInfo& rInfo = regions_[regionI];

            // Re-write mesh if dynamic or first time
            if (updateMesh)
            {
                meshWrite(rInfo);
            }

            // Write field data
            fieldWrite(rInfo);

            // Write cloud data
            cloudWrite(rInfo);
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
