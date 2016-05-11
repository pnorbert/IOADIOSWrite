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
#include "CStringList.H"

#include "FlatListOutput.H"
#include "ParticleBinaryBlob.H"

#include "adiosCloudHeaders.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::word Foam::adiosWrite::supportedCloudType(const word& cloudType)
{
    if (isNull(cloudType) || cloudType.empty())
    {
        return word::null;
    }

    // this needs reworking:
    if (cloudType == Cloud<indexedParticle>::typeName)
    {
        return Cloud<indexedParticle>::typeName;
    }
    else if (cloudType == Cloud<passiveParticle>::typeName)
    {
        return Cloud<passiveParticle>::typeName;
    }
    else if (cloudType == Cloud<basicKinematicCollidingParcel>::typeName)
    {
        // Not implemented: - use alternative
        return Cloud<basicKinematicParcel>::typeName;
    }
    else if (cloudType == Cloud<basicKinematicMPPICParcel>::typeName)
    {
        // Not implemented: - use alternative
        return Cloud<basicKinematicParcel>::typeName;
    }
    else if (cloudType == Cloud<basicKinematicParcel>::typeName)
    {
        return Cloud<basicKinematicParcel>::typeName;
    }

    // else if (cloudType == Cloud<basicReactingMultiphaseParcel>::typeName)
    // else if (cloudType == Cloud<basicReactingParcel>::typeName)
    // else if (cloudType == Cloud<basicThermoParcel>::typeName)
    // else if (cloudType == Cloud<molecule>::typeName)
    // else if (cloudType == Cloud<solidParticle>::typeName)
    // else if (cloudType == Cloud<basicSprayParcel>::typeName)

    return word::null;
}


size_t Foam::adiosWrite::cloudDefine(regionInfo& r)
{
    Info<< "  adiosWrite::cloudDefine: " << r.info() << endl;

    const fvMesh& mesh = time_.lookupObject<fvMesh>(r.name());

    r.cloudInfo_.clear();

    HashTable<const cloud*> allClouds = mesh.lookupClass<cloud>();

    HashTable<cloudInfo> cloudsUsed;
    HashTable<word> unsupported;

    forAllConstIter(HashTable<const cloud*>, allClouds, iter)
    {
        const word& name = (*iter)->name();
        if (findStrings(r.ignoredClouds_, name))
        {
            continue;
        }

        const word& type = (*iter)->type();

        Info<< "check cloud: " << name << " type = " << type << nl;

        if (r.autoWrite() || findStrings(r.requestedClouds_, name))
        {
            word useType = supportedCloudType(type);
            if (useType.empty())
            {
                unsupported.set(name, type);
            }
            else
            {
                cloudsUsed.insert
                (
                    name,
                    cloudInfo(name, type, useType)
                );
            }
        }
    }


    if (!unsupported.empty())
    {
        wordList names = unsupported.sortedToc();
        wordList types(names.size());

        forAll(names, nameI)
        {
            types[nameI] = unsupported[names[nameI]];
        }

        WarningInFunction
            << nl
            << unsupported.size() << " clouds not handled by adiosWrite" << nl
            << "  names: " << FlatListOutput<word>(names) << nl
            << "  types: " << FlatListOutput<word>(types) << nl << endl;
    }


    // Info<< "got clouds " << cloudsUsed << nl;

    size_t maxLen = 0;

    const wordList cloudNames = cloudsUsed.sortedToc();
    forAll(cloudNames, cloudI)
    {
        cloudInfo& cInfo = cloudsUsed[cloudNames[cloudI]];

        const word& cloudName = cInfo.name();
        const word& dispatch  = cInfo.dispatch();
        const fileName varPath = r.cloudPath(cInfo);

        regIOobject* obj = mesh.find(cloudName)();

        size_t bufLen = 0;

        // this needs reworking:
        if (dispatch == Cloud<indexedParticle>::typeName)
        {
            bufLen = cloudDefine
            (
                static_cast<Cloud<indexedParticle>&>(*obj),
                cInfo,
                varPath
            );
        }
        else if (dispatch == Cloud<passiveParticle>::typeName)
        {
            bufLen = cloudDefine
            (
                static_cast<Cloud<passiveParticle>&>(*obj),
                cInfo,
                varPath
            );
        }
        else if (dispatch == Cloud<basicKinematicCollidingParcel>::typeName)
        {
            bufLen = cloudDefine
            (
                static_cast<Cloud<basicKinematicCollidingParcel>&>(*obj),
                cInfo,
                varPath
            );
        }
        else if (dispatch == Cloud<basicKinematicMPPICParcel>::typeName)
        {
            bufLen = cloudDefine
            (
                static_cast<Cloud<basicKinematicMPPICParcel>&>(*obj),
                cInfo,
                varPath
            );
        }
        else if (dispatch == Cloud<basicKinematicParcel>::typeName)
        {
            bufLen = cloudDefine
            (
                static_cast<Cloud<basicKinematicParcel>&>(*obj),
                cInfo,
                varPath
            );
        }
        // else if (dispatch == Cloud<basicReactingMultiphaseParcel>::typeName)
        // else if (dispatch == Cloud<basicReactingParcel>::typeName)
        // else if (dispatch == Cloud<basicThermoParcel>::typeName)
        // else if (dispatch == Cloud<molecule>::typeName)
        // else if (dispatch == Cloud<solidParticle>::typeName)
        // else if (dispatch == Cloud<basicSprayParcel>::typeName)

        if (bufLen)
        {
            maxLen = Foam::max(maxLen, bufLen);
            r.cloudInfo_.append(cInfo);
        }
    }

    label nClouds = r.cloudInfo_.size();
    if (nClouds)
    {
        DynamicList<word> names(nClouds);
        forAllConstIter(SLList<cloudInfo>, r.cloudInfo_, iter)
        {
            const cloudInfo& cInfo = iter();

            names.append(cInfo.name());
        }

        const fileName varPath = r.regionPath();

        // number of active clouds as region attribute
        defineIntAttribute("nClouds", varPath, names.size());
        defineListAttribute("clouds", varPath, names);
    }

    return maxLen;
}


void Foam::adiosWrite::cloudWrite(const regionInfo& rInfo)
{
    const fvMesh& mesh = time_.lookupObject<fvMesh>(rInfo.name());
    Info<< "  adiosWrite::cloudWrite: " << rInfo.info() << endl;

    DynamicList<label> labelBuffer;
    DynamicList<scalar> scalarBuffer;

    forAllConstIter(SLList<cloudInfo>, rInfo.cloudInfo_, iter)
    {
        const cloudInfo& cInfo = iter();
        const word&  cloudName = cInfo.name();
        const word&   dispatch = cInfo.dispatch();

        // If the cloud contains no particles, jump to the next cloud
        // - this is probably redundant here
        if (cInfo.nTotal() == 0)
        {
            continue;
        }

        regIOobject* obj = mesh.find(cloudName)();

        // this needs reworking:
        if (dispatch == Cloud<indexedParticle>::typeName)
        {
            cloudWrite
            (
                static_cast<Cloud<indexedParticle>&>(*obj),
                cInfo,
                rInfo
            );
        }
        else if (dispatch == Cloud<passiveParticle>::typeName)
        {
            cloudWrite
            (
                static_cast<Cloud<passiveParticle>&>(*obj),
                cInfo,
                rInfo
            );
        }
        else if (dispatch == Cloud<basicKinematicCollidingParcel>::typeName)
        {
            cloudWrite
            (
                static_cast<Cloud<basicKinematicCollidingParcel>&>(*obj),
                cInfo,
                rInfo
            );
        }
        else if (dispatch == Cloud<basicKinematicMPPICParcel>::typeName)
        {
            cloudWrite
            (
                static_cast<Cloud<basicKinematicMPPICParcel>&>(*obj),
                cInfo,
                rInfo
            );
        }
        else if (dispatch == Cloud<basicKinematicParcel>::typeName)
        {
            cloudWrite
            (
                static_cast<Cloud<basicKinematicParcel>&>(*obj),
                cInfo,
                rInfo
            );
        }
        // else if (dispatch == Cloud<basicReactingMultiphaseParcel>::typeName)
        // else if (dispatch == Cloud<basicReactingParcel>::typeName)
        // else if (dispatch == Cloud<basicThermoParcel>::typeName)
        // else if (dispatch == Cloud<molecule>::typeName)
        // else if (dispatch == Cloud<solidParticle>::typeName)
        // else if (dispatch == Cloud<basicSprayParcel>::typeName)
    }

}


// ************************************************************************* //
