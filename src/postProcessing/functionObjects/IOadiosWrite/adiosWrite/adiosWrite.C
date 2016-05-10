/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 Norbert Podhorszki
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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
#include "HashSet.H"
#include "scalar.H"
#include "FlatListOutput.H"

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


void Foam::adiosWrite::regionInfo::read
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& topDict
)
{
    Info<< " Region dict: " << name_ << " (mesh " << mesh.name() << ")" << endl;

    if (!topDict.empty())
    {
        getAutoWrite(topDict); // check for auto-write
    }

    getAutoWrite(dict);  // check again for auto-write

    wordList toc(dict.toc());

    if (dict.found("objectNames") || !autoWrite_)
    {
        dict.lookup("objectNames") >> objectNames_;
    }
    if (dict.found("cloudNames") || !autoWrite_)
    {
        dict.lookup("cloudNames") >> cloudNames_;
    }
    if (dict.found("cloudAttribs") || !autoWrite_)
    {
        dict.lookup("cloudAttribs") >> cloudAttribs_;
    }

    // Do a basic check to see if the objectNames_ is accessible

    DynamicList<word> missing(objectNames_.size());
    forAll(objectNames_, i)
    {
        if (mesh.foundObject<regIOobject>(objectNames_[i]))
        {
            Info<< " " << objectNames_[i];
        }
        else
        {
            missing.append(objectNames_[i]);
        }
    }

    if (missing.size())
    {
        WarningInFunction
            << missing.size() << " objects not found in database:" << nl
            << "   " << FlatListOutput<word>(missing) << endl;
    }

    // Also print the cloud names
    forAll(cloudNames_, i)
    {
        Info<< ' ' << cloudNames_[i];
    }
    Info<< nl << endl;
}


size_t Foam::adiosWrite::defineVars(bool updateMesh)
{
    Info<< "adiosWrite::defineVars(" << (updateMesh ? "updateMesh" : "")
        << ") called at time "
        << obr_.time().timeName() << " time index "
        << obr_.time().timeIndex() << endl;

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
    defineIntAttribute("nRegions",  adiosCore::foamAttribute, regions_.size());

    // region-names attribute (list of strings)
    {
        stringList names(regions_.size());
        forAll(regions_, regionI)
        {
            names[regionI] = regions_[regionI].name_;
        }

        defineListAttribute("regions", adiosCore::foamAttribute, names);
    }

    // General information (as variable)
    defineIntVariable(adiosCore::timeAttribute/"index");
    defineScalarVariable(adiosCore::timeAttribute/"value");
    defineScalarVariable(adiosCore::timeAttribute/"deltaT");
    defineScalarVariable(adiosCore::timeAttribute/"deltaT0");

    forAll(regions_, regionI)
    {
        regionInfo& rInfo = regions_[regionI];
        const fileName varPath = rInfo.regionPath();

        const fvMesh& mesh = time_.lookupObject<fvMesh>(rInfo.name_);
        const polyPatchList& patches = mesh.boundaryMesh();

        stringList pNames(patches.size());
        stringList pTypes(patches.size());
        forAll(patches, patchI)
        {
            const polyPatch& p = patches[patchI];

            pNames[patchI] = p.name();
            pTypes[patchI] = p.type();
        }

        defineIntAttribute("nPatches",     varPath, patches.size());
        defineListAttribute("patch-names", varPath, pNames);
        defineListAttribute("patch-types", varPath, pTypes);

        bufLen = meshDefine(mesh, updateMesh);
        maxLen = Foam::max(maxLen, bufLen);

        bufLen = fieldDefine(rInfo);
        maxLen = Foam::max(maxLen, bufLen);

        bufLen = cloudDefine(rInfo);
        maxLen = Foam::max(maxLen, bufLen);
    }

    return maxLen;
}


void Foam::adiosWrite::deleteDefinitions()
{
    // Info<< "adiosWrite::deleteDefinitions() called at time "
    //     << obr_.time().timeName() << " time index "
    //     << obr_.time().timeIndex() << endl;

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
    obr_(obr),
    primaryMesh_(refCast<const fvMesh>(obr)),
    time_(primaryMesh_.time())
{
    Info<< "adiosWrite constructor called (" << Pstream::nProcs()
        << " procs) with primary mesh " << primaryMesh_.name() << endl;

    if (Pstream::nProcs() == 1)
    {
        // is serial - must do MPI_Init() ourselves
        MPI_Init(NULL, NULL);  // NULL args are OK for openmpi, mpich-2
    }

    MPI_Comm_dup(MPI_COMM_WORLD, &comm_);

    // Initialize ADIOS
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
        writeMethod_.c_str(),
        writeParams_.c_str(),
        ""  // base-path (unused) needs empty string, not a NULL pointer
    );

    // Write initial conditions (including mesh) if restartTime not set
    if (restartTime_ == VGREAT)
    {
        Info<< "write initial conditions" << endl;
        write();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::adiosWrite::~adiosWrite()
{
    adios_free_group(groupID_); // not necessary but nice to cleanup
    adios_finalize(Pstream::myProcNo());

    MPI_Comm_free(&comm_);

    if (Pstream::nProcs() == 1)
    {
        // is serial - must do MPI_Finalize() ourselves
        MPI_Finalize();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::adiosWrite::read(const dictionary& dict)
{
    // wordList toc_(dict.toc());
    // forAll(toc_, i)
    // {
    //     Info<< " TOC "<< i<<": " << toc_[i] << endl;
    // }
    // Info<< endl << endl;


    // test with auto-write
    Switch autoWrite = true;
    autoWrite.readIfPresent("autoWrite", dict); // top-level auto-write?

    // all known regions
    wordHashSet regionNames(time_.names<fvMesh>());
    // Info<< "Known regions: " << regionNames << endl;

    // Get the list of regions in adiosDict (as name list and dictionary list)
    const entry* entryPtr = dict.lookupEntryPtr("regions", false, false);

    label nRegion = 0;
    if (entryPtr->isDict())
    {
        // Regions dictionary
        const dictionary& allRegionsDict = entryPtr->dict();

        regions_.setSize
        (
            Foam::max
            (
                regionNames.size(),
                allRegionsDict.toc().size()
            )
        );

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
            const word& regName = iter().keyword();

            if (regionDict.lookupOrDefault<Switch>("active", true))
            {
                if (time_.foundObject<fvMesh>(regName))
                {
                    regions_[nRegion].name_  = regName;

                    regions_[nRegion].read
                    (
                        time_.lookupObject<fvMesh>(regName),
                        regionDict,
                        dict
                    );

                    Info<< regions_[nRegion].info() << endl;

                    ++nRegion;
                }
                else
                {
                    Info<<"no such region: " << regName << endl;
                }
            }

            // remove as treated
            regionNames.erase(regName);
        }
    }
    else if (autoWrite)
    {
        regions_.setSize(regionNames.size());
    }
    else if (!autoWrite)
    {
        FatalIOErrorInFunction(dict)
            << "regions must be specified in dictionary format"
            << exit(FatalIOError);
    }

    if (autoWrite)
    {
        // primary mesh first
        word regName = primaryMesh_.name();
        if (regionNames.found(regName))
        {
            regions_[nRegion].name_      = regName;
            regions_[nRegion].autoWrite_ = true;
            ++nRegion;

            // remove as treated
            regionNames.erase(regName);
        }

        wordList regNames = regionNames.sortedToc();
        forAll(regNames, regI)
        {
            const word& regName = regNames[regI];

            if (time_.foundObject<fvMesh>(regName))
            {
                // this lookup cannot actually fail
                regions_[nRegion].name_      = regName;
                regions_[nRegion].autoWrite_ = true;
                ++nRegion;
            }
        }
    }

    regions_.setSize(nRegion);

    // Default values prior to lookup in dictionary
    readMethod_  = "BP";
    writeMethod_ = "MPI";
    writeParams_ = "";

    dict.readIfPresent("readMethod",  readMethod_);
    dict.readIfPresent("adiosMethod", writeMethod_);  // old name
    dict.readIfPresent("writeMethod", writeMethod_);
    dict.readIfPresent("methodparams", writeParams_); // old name
    dict.readIfPresent("writeOptions", writeParams_);

    writeInterval_ = dict.lookupOrDefault<label>("writeInterval", 1);

    // Print info to terminal
    Info<< type() << " " << name() << ":" << endl
        << "  Compiled with " << adiosTraits<scalar>::adiosSize
        << " bytes precision." << nl
        << "  writing every " << writeInterval_ << " iterations:" << endl;


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

    Info<< type() << " " << name() << ":" << nl
        << "  ADIOS writeMethod: " << writeMethod_ << nl
        << "        writeParams: " << writeParams_ << nl
        << "     write interval: " << writeInterval_ << endl;
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

        // Info<<"times: " << times << endl;
        // Info<<"look for times in " << dataDirectory << endl;

        const instantList adiosTimes = adiosCore::findTimes();

        // Info<<"adios-times: " << adiosTimes << endl;
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
    Info<< "adiosWrite::end() called at time "
        << obr_.time().timeName() << " time index "
        << obr_.time().timeIndex() << endl;
}


void Foam::adiosWrite::timeSet()
{
    // Nothing to be done here
    Info<< "adiosWrite::timeSet() called at time "
        << obr_.time().timeName() << " time index "
        << obr_.time().timeIndex() << endl;
}


void Foam::adiosWrite::write()
{
    Info<< "adiosWrite::write() called at time "
        << obr_.time().timeName() << " time index "
        << obr_.time().timeIndex() << endl;

    size_t maxLen = 0;
    size_t bufLen = 0;

    // Check if we are going to write
    //if ( timeSteps_ == 0 )
    if (timeSteps_ == nextWrite_)
    {
        // Write info to terminal
        Info<< "Writing ADIOS data for time " << obr_.time().timeName() << endl;

        // Re-write mesh if first time or any of the meshes changed
        bool updateMesh = (timeSteps_ == 0);
        if (!updateMesh)
        {
            HashTable<const fvMesh*> allMeshes = time_.lookupClass<fvMesh>();

            forAllConstIter(HashTable<const fvMesh*>, allMeshes, iter)
            {
                updateMesh = (*iter)->changing();
                if (updateMesh)
                {
                    break;
                }
            }
        }

#if 0
        /* FIXME: code snippet to wait gdb attach on rank 0 process */
        /* if (timeSteps_ == 0)
        {
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
        }
        */
#endif

        // remove old ADIOS variables/attributes
        // - new variables may have appeared (eg, via function objects) etc
        // - patch sizes may change with every step
        // -> simply update everything for simplicity
        deleteDefinitions();

        classifyFields();

        // ADIOS requires to define all variables before writing anything
        Info<< "Define variables in ADIOS" << endl;
        bufLen = defineVars(updateMesh);
        maxLen = Foam::max(maxLen, bufLen);

        // Pout<<"reserve write buffer " << maxLen << endl;
        iobuffer_.reserve(maxLen); // NEEDS attention if iobuffer is not char!

        // Create/reopen ADIOS output file, and tell ADIOS how many bytes we
        // are going to write
        open();

        // General information (as variable)
        writeIntVariable
        (
            adiosCore::timeAttribute/"index",
            time_.timeIndex()
        );
        writeScalarVariable
        (
            adiosCore::timeAttribute/"value",
            time_.timeOutputValue()
        );
        writeScalarVariable
        (
            adiosCore::timeAttribute/"deltaT",
            time_.deltaT().value()
        );
        writeScalarVariable
        (
            adiosCore::timeAttribute/"deltaT0",
            time_.deltaT0().value()
        );

        forAll(regions_, regionI)
        {
            regionInfo& rInfo = regions_[regionI];

            const fvMesh& mesh = time_.lookupObject<fvMesh>(rInfo.name_);

            // Re-write mesh if dynamic or first time
            meshWrite(mesh, updateMesh);

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
