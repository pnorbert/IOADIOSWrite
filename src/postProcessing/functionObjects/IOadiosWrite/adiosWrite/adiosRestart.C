/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright(C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |              2015 Norbert Podhorszki
                            |              2016 OpenCFD Ltd.
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
#include "adiosReader.H"

#include "dictionary.H"
#include "scalar.H"
#include "basicKinematicCloud.H"
#include "emptyFvPatchField.H"

#include "IBufStream.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::adiosWrite::readData(const fileName& bpFile)
{
    Info<< " Read data of step " << bpFile << endl;

    adiosReader::helper helper(iobuffer_);

    if (!helper.open(bpFile.c_str(), comm_))
    {
        return false;
    }

    // Process N will read writeblock N
    helper.select(adios_selection_writeblock(Pstream::myProcNo()));
    helper.scan(true);
    helper.buffer.reserve(helper.maxLen);

    // direct lookup via time_.lookupClass<fvMesh>() would be nice,
    // but have to trick the compiler not to get the const-version.
    // - instead, get the names and do the lookup ourselves

    wordList meshNames = time_.names<fvMesh>();
    sort(meshNames);

    // this is somewhat like demand-driven loading
    forAll(meshNames, regI)
    {
        const word& regName = meshNames[regI];
        fvMesh& mesh = const_cast<fvMesh&>(time_.lookupObject<fvMesh>(regName));

        Info<<"lookup: " << regName << endl;

        HashTable<adiosReader::fieldInfo> fromFile = helper.getFieldInfo(regName);
        wordList fieldNames = fromFile.sortedToc();

        forAll(fieldNames, fieldI)
        {
            const word& name = fieldNames[fieldI];
            const adiosReader::fieldInfo& src = fromFile[name];

            if ((static_cast<objectRegistry&>(mesh)).found(name))
            {
                readVolField(mesh.find(name)(), helper, src);
            }
        }
    }


    bool ok = true;
    forAll(regions_, regionId)
    {
        regionInfo& rInfo = regions_[regionId];

        ok = readClouds(helper, rInfo);
    }

    helper.close();
    return ok;
}


bool Foam::adiosWrite::readData(const instant& when)
{
    Info<< " Read data of step " << when.name() << endl;
    return readData(when.name());
}


bool Foam::adiosWrite::readData()
{
    Info<< " Read data of step " << obr_.time().timeName() << endl;
    return readData(dataDirectory/obr_.time().timeName() + ".bp");
}



bool Foam::adiosWrite::readVolField
(
    regIOobject* obj,
    adiosReader::helper& helper,
    const adiosReader::fieldInfo& src
)
{
    const word& fieldType = obj->type();
    const word&   srcType = src.type();

    if (fieldType != srcType)
    {
        // probably fatal:
        Info<<"WARNING mismatch on field " << src
            << " (expected: " << fieldType
            << " but had " << srcType << ")\n";
        return false;
    }

    Info<<"read " << obj->name() << " (type " << fieldType << ") from file\n";
    if (fieldType == volScalarField::typeName)
    {
        return fieldRead
        (
            static_cast<volScalarField&>(*obj),
            helper,
            src
        );
    }
    else if (fieldType == volVectorField::typeName)
    {
        return fieldRead
        (
            static_cast<volVectorField&>(*obj),
            helper,
            src
        );
    }
    else if (fieldType == surfaceScalarField::typeName)
    {
        return fieldRead
        (
            static_cast<surfaceScalarField&>(*obj),
            helper,
            src
        );
    }
    else if (fieldType == volSphericalTensorField::typeName)
    {
        return fieldRead
        (
            static_cast<volSphericalTensorField&>(*obj),
            helper,
            src
        );
    }
    else if (fieldType == volSymmTensorField::typeName)
    {
        return fieldRead
        (
            static_cast<volSymmTensorField&>(*obj),
            helper,
            src
        );
    }
    else if (fieldType == volTensorField::typeName)
    {
        return fieldRead
        (
            static_cast<volTensorField&>(*obj),
            helper,
            src
        );
    }
    else
    {
        return false;
    }
}


bool Foam::adiosWrite::readClouds(adiosReader::helper& helper, regionInfo& rInfo)
{
    bool ok = true;
    const fvMesh& mesh = time_.lookupObject<fvMesh>(rInfo.name_);

    forAll(rInfo.cloudNames_, cloudI)
    {
        Info<< "    cloud: " << rInfo.cloudNames_[cloudI] << endl;

        const kinematicCloud& constCloud =
            mesh.lookupObject<kinematicCloud>(rInfo.cloudNames_[cloudI]);

        kinematicCloud& cloud = const_cast<kinematicCloud&>(constCloud);

        //basicKinematicCloud *q =(basicKinematicCloud*) &cloud;
        basicKinematicCloud *q = reinterpret_cast<basicKinematicCloud*>(&cloud);

        fileName varPath = rInfo.cloudPath(rInfo.cloudNames_[cloudI]);
        fileName datasetName(varPath/"nParticle");

        // Get the number of particles saved by this rank
        int nparts = 0;
        ok = helper.getDataSet(datasetName, &nparts);

        /*
        // Get the number of particles saved by this rank
        ADIOS_VARINFO * vi = adios_inq_var(helper.file, datasetName);
        if (vi != NULL)
        {
            adios_inq_var_blockinfo(helper.file, vi);
            if (vi->sum_nblocks > Pstream::myProcNo())
            {
                nParticles = vi->blockinfo[Pstream::myProcNo()].count[0];
            }
            else
            {
                WarningInFunction
                    << "Error reading cloud " << datasetName
                    << ": Number of available blocks = " << vi->sum_nblocks
                    << " is less than this rank " << Pstream::myProcNo()
                    << endl;
            }
            adios_free_varinfo(vi);
        }
        else
        {
            WarningInFunction
                << "Error reading cloud " << datasetName
                << " from adios checkpoint file: variable not found."
                << endl;
            ok = false;
        }
        */

        if (!ok) break;

////        // Read into a plain continuous array for the data
////        // Allocate memory for 1-comp. dataset of type 'integer' for adios_integer reads
////        int labelData[nparts];
////
////        // Read original processor ID
////        if (findStrings(rInfo.cloudAttribs_, "origProc"))
////        {
////            Info<< "      dataset origProc " << endl;
////            ok = helper.getDataSet(varPath/"origProc", labelData);
////            label i = 0;
////            forAllIter(basicKinematicCloud, *q, pIter)
////            {
////                pIter().origProc() = labelData[i++];
////            }
////        }
////
        if (!ok) break;
    }

    return ok;
}


// ************************************************************************* //
