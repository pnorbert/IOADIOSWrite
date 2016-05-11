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
#include "hashedWordList.H"
#include "HashSet.H"
#include "scalar.H"
#include "SortableList.H"
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //


namespace Foam
{
//! \cond fileScope
static dictionary subOrEmptyDict
(
    const dictionary& dict,
    const word& keyword
)
{
    const entry* entryPtr = dict.lookupEntryPtr(keyword, false, false);
    if (entryPtr)
    {
        // let it fail if it is the incorrect type
        return entryPtr->dict();
    }
    else
    {
        return dictionary(dict, dictionary(dict.name() + '.' + keyword));
    }
}
}
//! \endcond


// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

void Foam::adiosWrite::regionInfo::read
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    resetAll();
    name_ = mesh.name(); // ensure absolute consistency

    Info<< " Region dict: " << name_ << endl;

    const dictionary* writePtr  = dict.subDictPtr("write");
    const dictionary* ignorePtr = dict.subDictPtr("ignore");

    const dictionary& wrtDict = (writePtr  ? *writePtr  : dictionary::null);
    const dictionary& ignDict = (ignorePtr ? *ignorePtr : dictionary::null);

    explicitWrite_.readIfPresent("explicit", wrtDict);

    if (explicitWrite_ || wrtDict.found("fields"))
    {
        wrtDict.lookup("fields") >> requestedFields_;
    }
    if (explicitWrite_ || wrtDict.found("clouds"))
    {
        wrtDict.lookup("clouds") >> requestedClouds_;
    }
    if (explicitWrite_ || wrtDict.found("cloudAttrs"))
    {
        wrtDict.lookup("cloudAttrs") >> requestedAttrs_;
    }

    if (ignDict.found("fields"))
    {
        ignDict.lookup("fields") >> ignoredFields_;
    }
    if (ignDict.found("clouds"))
    {
        ignDict.lookup("clouds") >> ignoredClouds_;
    }

    // Check if the requested fields are actually accessible
    DynamicList<word> missing(requestedFields_.size());
    forAll(requestedFields_, i)
    {
        const wordRe& what = requestedFields_[i];

        if (!what.isPattern())
        {
            if (mesh.foundObject<regIOobject>(what))
            {
                Info<< " " << what;
            }
            else
            {
                missing.append(what);
            }
        }
    }

    // Also print the cloud names
    forAll(requestedClouds_, i)
    {
        Info<< ' ' << requestedClouds_[i];
    }

    Info<< nl << endl;

    if (missing.size())
    {
        WarningInFunction
            << nl
            << missing.size() << " objects not found in database:" << nl
            << "   " << FlatListOutput<word>(missing) << endl;
    }
}


void Foam::adiosWrite::regionInfo::resetAll()
{
    explicitWrite_ = false;
    requestedFields_.clear();
    requestedClouds_.clear();
    requestedAttrs_.clear();
    ignoredFields_.clear();
    ignoredClouds_.clear();
    fieldsToWrite_.clear();
    cloudInfo_.clear();
}


size_t Foam::adiosWrite::defineVars()
{
    Info<< "adiosWrite::defineVars() called at time "
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

    // other general information
    defineIntAttribute("nProcs",    adiosCore::foamAttribute, Pstream::nProcs());
    defineIntAttribute("nRegions",  adiosCore::foamAttribute, regions_.size());

    // region-names attribute (list of strings)

    // mesh updated - may also want moved points etc.
    // this is still not the greatest:
    {
        DynamicList<word> names(regions_.size());
        DynamicList<word> topo(regions_.size());
        DynamicList<word> moved(regions_.size());

        forAllConstIter(SLList<regionInfo>, regions_, iter)
        {
            const regionInfo& rInfo = iter();
            const word& regName = rInfo.name();

            names.append(regName);

            if (rInfo.topoChanging())
            {
                topo.append(rInfo.name());
            }
            else if (rInfo.moving())
            {
                moved.append(rInfo.name());
            }
        }

        defineListAttribute("regions",     adiosCore::foamAttribute, names);
        defineListAttribute("moving",      adiosCore::foamAttribute, moved);
        defineListAttribute("topo-change", adiosCore::foamAttribute, topo);
    }

    // General information (as variable)
    defineIntVariable(adiosTime::attr[adiosTime::INDEX]);
    defineScalarVariable(adiosTime::attr[adiosTime::VALUE]);
    defineScalarVariable(adiosTime::attr[adiosTime::DT]);
    defineScalarVariable(adiosTime::attr[adiosTime::DT0]);

    forAllIter(SLList<regionInfo>, regions_, iter)
    {
        regionInfo& rInfo = iter();
        const fileName varPath = rInfo.regionPath();

        const fvMesh& mesh = time_.lookupObject<fvMesh>(rInfo.name());
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

        bufLen = meshDefine(rInfo);
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
    time_(refCast<const fvMesh>(obr).time()),
    restartIndex_(-1)
{
    Info<< "adiosWrite constructor called (" << Pstream::nProcs()
        << " procs) with primary mesh "
        << refCast<const fvMesh>(obr).name() << endl;

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


Foam::adiosWrite::regionInfo::regionInfo
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    name_(mesh.name()),
    explicitWrite_(false),
    requestedFields_(),
    ignoredFields_(),
    topoState_(polyMesh::UNCHANGED),
    requestedClouds_(),
    ignoredClouds_(),
    requestedAttrs_(),
    fieldsToWrite_(),
    cloudInfo_()
{
    read(mesh, dict);
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
    regions_.clear();

    // all known regions
    // - in sorted order for consistency, hashed for quick lookup
    hashedWordList regionNames;
    {
        SortableList<word> sorted(time_.names<fvMesh>());
        regionNames.transfer(sorted);
    }

    // Info<< "Known regions: " << regionNames << endl;

    // top-level dictionaries
    dictionary regionsDict = subOrEmptyDict(dict, "regions");

    // verify that "regions" only contains sub-dictionaries
    // and detect any unknown region names
    {
        wordList bad(regionsDict.size());
        label nbad = 0;

        wordList warn(regionsDict.size());
        label nwarn = 0;

        forAllConstIter(dictionary, regionsDict, iter)
        {
            const word& regName = iter().keyword();

            if (!regionNames.contains(regName))
            {
                warn[nwarn++] = regName;
            }
            if (!iter().isDict())
            {
                bad[nbad++] = regName;
            }
        }

        if (nwarn)
        {
            warn.setSize(nwarn);

            WarningInFunction
                << warn.size() << " unknown regions specified in adiosDict:"
                << nl
                << "   " << FlatListOutput<word>(warn) << endl;
        }

        if (nbad)
        {
            bad.setSize(nbad);

            FatalErrorInFunction
                << bad.size() << " bad region specifications in adiosDict:"
                << nl
                << "   " << FlatListOutput<word>(bad) << nl
                << "must be dictionary format with region-name as the keyword"
                << nl
                << exit(FatalIOError);
        }
    }

    // for each region
    // - uses "explicit" entry and "write" and "ignore" sub-dictionaries
    forAll(regionNames, regI)
    {
        // cannot fail since the region name came from previous lookup
        const word& regName = regionNames[regI];
        const fvMesh& mesh = time_.lookupObject<fvMesh>(regName);

        const entry* dictPtr = regionsDict.lookupEntryPtr(regName, false, false);

        if (dictPtr)
        {
            // make a copy and merge in region-specific content
            // (overwrite existing entries)

            dictionary mergedDict(regionsDict);
            mergedDict <<= dictPtr->dict();

            regions_.append(regionInfo(mesh, mergedDict));
        }
        else
        {
            regions_.append(regionInfo(mesh, dict));
        }

        Info<< regions_.last().info() << endl;
    }


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


bool Foam::adiosWrite::restart()
{
    static bool restarted = false;

    if (restarted || restartTime_ == VGREAT)
    {
        return false;
    }
    restarted = true; // even if it fails don't try again

    Info<< "adiosWrite::restart(" << restartTime_ << ")" << endl;

    const instantList adiosTimes = adiosCore::findTimes();
    const label adiosTimeIndex = Time::findClosestTimeIndex(adiosTimes, restartTime_);

    if (adiosTimeIndex < 0)
    {
        FatalErrorInFunction
            << "No appropriate ADIOS restart found for time " << restartTime_
            << exit(FatalIOError);
    }

    Info<<"restart from adios " << adiosTimes[adiosTimeIndex] << endl;

    // Classify fields in the object space,
    // then read data from restart file for each
    classifyFields();

    adiosTime bpTime = readData(adiosTimes[adiosTimeIndex]);
    if (bpTime.valid())
    {
        Time& t = const_cast<Time&>(obr_.time());
        t.setTime(bpTime.timeValue(), bpTime.timeIndex());

        Info<< "Time = " << t.timeOutputValue() << " index "
            << t.timeIndex() << endl;

        restartIndex_ = bpTime.timeIndex();
        // TODO: handle deltaT, deltaT0 etc
    }
    else
    {
        FatalErrorInFunction
            << "Restart reading failed for time " << restartTime_
            << exit(FatalIOError);
    }

    return true; // even if it fails we don't try it again
}


void Foam::adiosWrite::execute()
{
    // execute() is called at the first non-zero time, after the calculation,
    // before write() is called.
    // This is a good point to do restart because fields created by other
    // function objects exist at this point (e.g. fieldAverage variables)

    Info<< "adiosWrite::execute() called at time "
        << obr_.time().timeName() << " time index "
        << obr_.time().timeIndex() << endl;

    static bool restarted = false;
    if (!restarted)
    {
        restarted = true;
        restart();
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

    const bool writeNow = (timeSteps_ == nextWrite_);
    if (writeNow)
    {
        // time of next write
        nextWrite_ = timeSteps_ + writeInterval_;
    }

    if (writeNow && restartIndex_ == obr_.time().timeIndex())
    {
        Info<< "Writing ADIOS data for time " << obr_.time().timeName()
            << " ... skipped at restart time" << endl;
    }
    else if (writeNow)
    {
        Info<< "Writing ADIOS data for time " << obr_.time().timeName() << endl;

        // Re-write mesh if first time or any of the meshes changed
        forAllConstIter(SLList<regionInfo>, regions_, iter)
        {
            const regionInfo& rInfo = iter();

            rInfo.meshChanging
            (
                (timeSteps_ == 0),
                time_.lookupObject<fvMesh>(rInfo.name())
            );
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

        classifyFields(true); // verbose

        // ADIOS requires to define all variables before writing anything
        Info<< "Define variables in ADIOS" << endl;

        size_t maxLen = defineVars();

        // Pout<<"reserve write buffer " << maxLen << endl;
        iobuffer_.reserve(maxLen); // NEEDS attention if iobuffer is not char!

        // Create/reopen ADIOS output file, and tell ADIOS how many bytes we
        // are going to write
        open();

        // General information (as variable)
        writeIntVariable
        (
            adiosTime::attr[adiosTime::INDEX],
            time_.timeIndex()
        );
        writeScalarVariable
        (
            adiosTime::attr[adiosTime::VALUE],
            time_.timeOutputValue()
        );
        writeScalarVariable
        (
            adiosTime::attr[adiosTime::DT],
            time_.deltaTValue()
        );
        writeScalarVariable
        (
            adiosTime::attr[adiosTime::DT0],
            time_.deltaT0Value()
        );

        forAllIter(SLList<regionInfo>, regions_, iter)
        {
            regionInfo& rInfo = iter();

            // Re-write mesh if dynamic or first time
            meshWrite(rInfo);

            // Write field data
            fieldWrite(rInfo);

            // Write cloud data
            cloudWrite(rInfo);
        }

        // Close the ADIOS dataset at every timestep. It must be closed!
        // Data is flushed from memory at this point
        close();
    }

    // Update counter
    ++timeSteps_;
}


// ************************************************************************* //
