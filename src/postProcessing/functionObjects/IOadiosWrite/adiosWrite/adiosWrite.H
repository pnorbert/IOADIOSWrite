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

Class
    Foam::adiosWrite

Description
    Writes fields and particle clouds to an ADIOS file.
    Based on the IOh5Write module created by
       Håkon Strandenes, hakostra@stud.ntnu.no

SourceFiles
    ParticleBinaryBlob.C
    ParticleBinaryBlobFragment.C
    adiosCloudTemplates.C
    adiosCore.C
    adiosFieldTemplates.C
    adiosRestart.C
    adiosWrite.C
    adiosWriteCloud.C
    adiosWriteField.C
    adiosWriteFunctionObject.C

Author
    Norbert Podhorszki, pnorbert@ornl.gov, 2015

\*---------------------------------------------------------------------------*/

#ifndef adiosWrite_H
#define adiosWrite_H

#include "adiosCore.H"
#include "adiosReader.H"
#include "adiosCoreWrite.H"
#include "adiosTime.H"

#include "runTimeSelectionTables.H"
#include "interpolation.H"
#include "fvMesh.H"

#include "adios.h"
#include "adios_read.h"

/**
 * expand the particle information into separate fields
 * useful for debugging
 */
// #define FOAM_ADIOS_CLOUD_EXPAND
// #undef  FOAM_ADIOS_CLOUD_EXPAND

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;

/*---------------------------------------------------------------------------*\
                         Class adiosWrite Declaration
\*---------------------------------------------------------------------------*/

class adiosWrite
:
    public adiosCoreWrite
{
    // Forward declarations of various sub-components
        class regionInfo;
        class cloudInfo;


    // Private data

        // Reference to the object registry (fvMesh)
        const objectRegistry& obr_;

        //- Reference to the time database
        const Time& time_;

        //- track time index from restart (avoid write immediately after restart)
        int restartIndex_;

        //- Restart time requested
        scalar restartTime_;

        //- List of region information
        typedef SLList<regionInfo> RegionInfoContainer;

        //- List of regions holding all data needed (mesh, names)
        RegionInfoContainer regions_;

        //- Write interval (in number of iterations)
        label writeInterval_;

        // Used to track when to write
        label nextWrite_;
        label timeSteps_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        adiosWrite(const adiosWrite&);

        //- Disallow default bitwise assignment
        void operator=(const adiosWrite&);

protected:

    // Static Member Functions

        //- Return name of cloud-type (or surrogate) or empty word
        static word supportedCloudType(const word& cloudType);


    // Protected Member Functions

        // General functions

            //- Create and open dataset for all data at given time-step
            bool open(const Time&);


        // Functions for handling field data

            //- Define fields for the region
            //  Return the number of fields defined.
            label fieldDefine(const regionInfo&);

            //- Field write wrapper
            void fieldWrite(const regionInfo&);

            //- Read data for a given volume field type
            bool readVolField
            (
                regIOobject* obj,
                const adiosReader& reader,
                const adiosReader::fieldInfo& src
            );

            //- Read data for a given volume field type
            bool readCloud
            (
                const fvMesh& mesh,
                const word& cloudName,
                const adiosReader& reader,
                const adiosReader::cloudInfo& src
            );

            //- Read data for a given type
            template<class FieldType>
            bool fieldRead
            (
                FieldType& field,
                const adiosReader& reader,
                const adiosReader::fieldInfo& src
            );

            //- Read data for a given type
            template<class FieldType>
            bool fieldRead
            (
                const adiosReader& reader,
                const fvMesh& mesh,
                const adiosReader::fieldInfo& src
            );


      // Functions for handling particle cloud data

            //- Define ADIOS cloud variables/attributes for the region.
            //  Return the number of clouds defined.
            label cloudDefine(regionInfo&);


            //- Cloud define for a given type
            template<class CloudType>
            bool cloudDefine
            (
                const CloudType& cloud,
                cloudInfo& cInfo,
                const fileName& varPath
            );


            //- Write ADIOS cloud variables
            void cloudWrite(const regionInfo&);

            //- Write ADIOS cloud variable for a given type
            template<class CloudType>
            void cloudWrite
            (
                const CloudType& cloud,
                const cloudInfo& cInfo,
                const regionInfo& rInfo
            );


        // Miscallenous definitions

            //- Define ADIOS variables (before writing them)
            //  Mesh information is defined if the region has been updated
            void defineVars();

            //- Write all known ADIOS variables
            void writeVars();


        // Restart functions

            adiosTime readData(const fileName&);
            adiosTime readData(const instant& when);

            //- Read data from current time
            //  Attempt to read data for all variables existing in memory
            //  Returns time information on success, invalid on failure
            adiosTime readData();

            //- Restart functionality
            bool restart();

public:

    //- Runtime type information
    TypeName("adiosWrite");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        adiosWrite
        (
            const word& groupName,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~adiosWrite();


    // Member Functions

        //- Read the adiosWrite data
        virtual void read(const dictionary&);

        //- Execute, handles restart at first call if requested
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++,
        //  Note: currently does nothing
        virtual void timeSet();

        //- Write a timestep to file
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&)
        {}
};


//- Per-region info variables grouped together to be able to create a list of them
class adiosWrite::cloudInfo
{
    //- Name of the cloud
    word name_;

    //- Type of the cloud
    word type_;

    //- Alternative dispatch type
    word dispatch_;

    //- Total number of particles (sum of nParticle_ list)
    label nTotal_;

    //- Array containing number of particles per process (often used list)
    List<label> nParticle_;

public:

    //- Construct with name and type
    cloudInfo(const word& name, const word& type)
    :
        name_(name),
        type_(type),
        dispatch_(),
        nTotal_(0),
        nParticle_(Pstream::nProcs())
    {
        reset();
    }


    cloudInfo(const word& name, const word& type, const word& useType)
    :
        name_(name),
        type_(type),
        dispatch_(),
        nTotal_(0),
        nParticle_(Pstream::nProcs())
    {
        reset();

        if (type != useType)
        {
            dispatch_ = useType;
        }
    }


    ~cloudInfo() {}


    //- Name of the cloud
    const word& name() const
    {
        return name_;
    }


    //- Type of the cloud
    const word& type() const
    {
        return type_;
    }


    //- Dispatch type for the cloud
    const word& dispatch() const
    {
        return dispatch_.empty() ? type_ : dispatch_;
    }


    //- Uses indirect support
    bool needsDispatch() const
    {
        return !dispatch_.empty();
    }


    //- Total number of particles across all processes
    label nTotal() const
    {
        return nTotal_;
    }


    //- The number of particles on this process
    label nParticle() const
    {
        return nParticle_[Pstream::myProcNo()];
    }


    //- Set the number of particles and sync across processes
    //  Return the overall total
    label nParticle(label n)
    {
        // The number of particles on each process
        nParticle_[Pstream::myProcNo()] = n;
        Pstream::gatherList(nParticle_);
        Pstream::scatterList(nParticle_);

        // Sum total number of particles on all processes
        nTotal_ = sum(nParticle_);

        return nTotal_;
    }


    //- Local process offset within a global array of parcels
    label offset() const
    {
        label off = 0;
        for (label prev=0; prev < Pstream::myProcNo(); ++prev)
        {
            off += nParticle_[prev];
        }

        return off;
    }


    //- Reset to zeros
    void reset()
    {
        nTotal_ = 0;
        forAll(nParticle_, i)
        {
            nParticle_[i] = 0;
        }
    }

};


//- Per-region info variables grouped together to be able to create a list of them
class adiosWrite::regionInfo
{
    //- Name of the region
    word name_;

    //- Normally don't want an explicit-write (want auto-write)
    Switch explicitWrite_;

    //- Names of volume fields to write
    wordReList requestedFields_;

    //- Names of volume fields to ignore
    wordReList ignoredFields_;

    // on-the-fly information
    mutable enum polyMesh::readUpdateState topoState_;


    //- Read info for one region from its own dictionary
    void read(const fvMesh& mesh, const dictionary& dict);

    //- reset all information
    void resetAll();


public:

    //- Names of clouds to write
    wordReList requestedClouds_;

    //- Names of clouds to ignored
    wordReList ignoredClouds_;

    //- Additional attributes to write for particle data
    wordReList requestedAttrs_;

    //- Fields to be written. Field name as key, class as value
    HashTable<word> fieldsToWrite_;

    //- List of cloud information
    typedef SLList<cloudInfo> CloudInfoContainer;

    //- Clouds to be processed
    CloudInfoContainer cloudInfo_;


    //- Construct null
    regionInfo()
    :
        name_(),
        explicitWrite_(false)
    {};


    //- Construct from dictionary information
    regionInfo
    (
        const fvMesh& mesh,
        const dictionary& dict
    );


    ~regionInfo() {};


    //- If auto-write is in effect
    const word& name() const
    {
        return name_;
    }


    //- If auto-write is in effect
    bool autoWrite() const
    {
        return !explicitWrite_;
    }


    //- Set topology changing information, optionally force update
    void meshChanging(const fvMesh& mesh, bool force = false) const
    {
        if (force)
        {
            topoState_ = polyMesh::TOPO_CHANGE;
        }
        else if (mesh.topoChanging())
        {
            topoState_ = polyMesh::TOPO_CHANGE;
        }
        else if (mesh.moving())
        {
            topoState_ = polyMesh::POINTS_MOVED;
        }
        else
        {
            topoState_ = polyMesh::UNCHANGED;
        }
    }


    //- Is mesh moving
    bool moving() const
    {
        return topoState_ == polyMesh::POINTS_MOVED;
    }


    //- Is mesh topology changing
    bool topoChanging() const
    {
        return
        (
            topoState_ == polyMesh::TOPO_CHANGE
         || topoState_ == polyMesh::TOPO_PATCH_CHANGE
        );
    }


    //- Is mesh changing (topology changing and/or moving)
    bool changing() const
    {
        return moving() || topoChanging();
    }


    //- Classify fields into the appropriate groups,
    //  returns the number of fields
    label classifyFields
    (
        const fvMesh& mesh,
        bool verbose = false
    );


    //- Path name for regions
    inline string regionPath() const
    {
        return adiosCore::regionPath(name_);
    }

    //- Path name for fields
    inline fileName fieldPath(const string& var = string::null) const
    {
        return adiosCore::fieldPath(name_, var);
    }


    //- Path name for named cloud
    inline fileName cloudPath(const string& cloudName) const
    {
        return adiosCore::cloudPath(name_, cloudName);
    }


    //- Path name for named cloud
    inline fileName cloudPath(const cloudInfo& cinfo) const
    {
        return cloudPath(cinfo.name());
    }


    //- Simple text representation for output purposes
    inline word info() const
    {
        return name_;
    }

};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "adiosCloudTemplates.C"
    #include "adiosFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
