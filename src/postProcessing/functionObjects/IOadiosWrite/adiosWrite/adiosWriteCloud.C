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

#include "ParticleBinaryBlob.H"
#include "basicKinematicCloud.H"
#include "reactingCloud.H"

#include <stdio.h> // sprintf

// * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * * //

// file-local
static inline Foam::fileName oldCloudPath(Foam::label index, std::string name)
{
    return "region" + Foam::name(index) / "clouds" / name;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

size_t Foam::adiosWrite::cloudDefine(regionInfo& r)
{
    size_t bufLen = 0;
    size_t maxLen = 0;

    Info<< "  adiosWrite::cloudDefine: region "
        << r.index_ << " " << r.name_ << ": " << endl;

    const fvMesh& mesh = time_.lookupObject<fvMesh>(r.name_);

#ifdef FOAM_ADIOS_CLOUD_EXPLICIT_NAMES
    cloudDefineExplicit(r);
#endif

    // HashTable<const cloud*> allClouds = mesh.lookupClass<cloud>();

    // Info<< cloud.type() << endl;
    // Info<<"clouds: " << allClouds << endl;

//     forAllConstIter(HashTable<const cloud*>, allClouds, iter)
//     {
//         Info<<"  cloud=" << (*iter)->name() << endl;
//         Info<<"  class=" << (*iter)->type() << endl;
//
//         if (isA<reactingCloud>(*(*iter)))
//         {
//             Info<<"  reacting-cloud YES\n";
//         }
//         else
//         {
//             Info<<"  reacting-cloud NO\n";
//         }
//     }

    label nClouds = 0;
    forAll(r.cloudNames_, cloudI)
    {
        Info<< "    cloud: " << r.cloudNames_[cloudI] << endl;

        const kinematicCloud& cloud =
            mesh.lookupObject<kinematicCloud>(r.cloudNames_[cloudI]);

        const basicKinematicCloud& q =
            static_cast<const basicKinematicCloud&>(cloud);

        // Number of particles on this process
        const label myParticles = q.size();

        // Find the number of particles on each process
        r.nParticles_[Pstream::myProcNo()] = myParticles;
        Pstream::gatherList(r.nParticles_);
        Pstream::scatterList(r.nParticles_);

        // Sum total number of particles on all processes
        r.nTotalParticles_ = sum(r.nParticles_);

        // If the cloud contains no particles, jump to the next cloud
        if (r.nTotalParticles_ == 0)
        {
            Info<< "    " << r.cloudNames_[cloudI]
                <<": No particles in cloud. Skipping definition."
                << endl;

            continue;
        }

        const fileName varPath = r.cloudPath(nClouds++);

        defineAttribute("class", varPath, cloud.type());  // cloud type
        defineAttribute("name",  varPath, q.name());      // cloud name

        // total number of particles as an attribute
        defineIntAttribute("nParticle", varPath, r.nTotalParticles_);

        // number of particles per processor as a field
        defineIntVariable(varPath/"nParticle");

        // determine blob sizes

        List<label> blobSize(Pstream::nProcs(), 0);
        if (myParticles)
        {
            ORawCountStream os;  // always raw binary content
            os << *(q.first());
            blobSize[Pstream::myProcNo()] = os.size();
        }
        Pstream::gatherList(blobSize);
        Pstream::scatterList(blobSize);

        Info<< "particles: " << r.nParticles_ << endl;
        Info<< "blob-size: " << blobSize << endl;

        ParticleBinaryBlob blob
        (
            basicKinematicCloud::particleType::propertyTypes(),
            basicKinematicCloud::particleType::propertyList()
        );

        Info<< "Blob Mapping\n" << blob << nl
            << "expected blob size " <<  blob.byteSize() << nl
            << blob.types() << nl << blob.names() << endl;

        // common - local offset value
        label localOffset = 0;
        for (label procI=1; procI < Pstream::nProcs(); ++procI)
        {
            localOffset += r.nParticles_[procI-1];
        }
        const string offsetDims = Foam::name(localOffset);

        // for scalars
        const string localDims  = Foam::name(myParticles);
        const string globalDims = Foam::name(r.nTotalParticles_);

        // for vectors
        const string localDims3  = localDims  + ",3";
        const string globalDims3 = globalDims + ",3";

        {
            // stream contents
            const fileName varName = varPath/"__blob__";

            const string localBlob  = localDims  + "," + Foam::name(blob.byteSize());
            const string globalBlob = globalDims + "," + Foam::name(blob.byteSize());

            adios_define_var
            (
                groupID_,
                varName.c_str(),                // name
                NULL,                           // path (deprecated)
                adios_unsigned_byte,            // data-type
                localBlob.c_str(),              // local dimensions
                globalBlob.c_str(),             // global dimensions
                offsetDims.c_str()              // local offsets
            );

            bufLen = myParticles * blob.byteSize();
            maxLen = Foam::max(maxLen, bufLen);

            outputSize_ += bufLen;

            defineListAttribute("names",     varName, blob.names());
            defineListAttribute("types",     varName, blob.types());
            defineListAttribute("offset",    varName, blob.offsets());
            defineListAttribute("byte-size", varName, blob.byteSizes());
        }

        // additional output
        label nLabels  = 0;
        label nScalars = 0;
        label nVectors = 0;

        // additional output - declare in case it is used
        {
            adios_define_var
            (
                groupID_,
                (varPath/"Us").c_str(),
                NULL,
                adiosTraits<scalar>::adiosType,
                localDims3.c_str(),
                globalDims3.c_str(),
                offsetDims.c_str()
            );
            ++nVectors;
        }


#ifdef FOAM_ADIOS_CLOUD_EXPAND
        // walk the blob framents
        {
            forAllConstIter(SLList<ParticleBinaryBlob::Fragment>, blob, iter)
            {
                const ParticleBinaryBlob::Fragment& frag = *iter;

                if (frag.type() == pTraits<label>::typeName)
                {
                    adios_define_var
                    (
                        groupID_,
                        (varPath/frag.name()).c_str(),
                        NULL,
                        adiosTraits<label>::adiosType,
                        localDims.c_str(),
                        globalDims.c_str(),
                        offsetDims.c_str()
                    );
                    ++nLabels;
                }
                else if (frag.type() == pTraits<scalar>::typeName)
                {
                    adios_define_var
                    (
                        groupID_,
                        (varPath/frag.name()).c_str(),
                        NULL,
                        adiosTraits<scalar>::adiosType,
                        localDims.c_str(),
                        globalDims.c_str(),
                        offsetDims.c_str()
                    );
                    ++nScalars;
                }
                else if (frag.type() == pTraits<vector>::typeName)
                {
                    adios_define_var
                    (
                        groupID_,
                        (varPath/frag.name()).c_str(),
                        NULL,
                        adiosTraits<scalar>::adiosType,
                        localDims3.c_str(),
                        globalDims3.c_str(),
                        offsetDims.c_str()
                    );
                    ++nVectors;
                }
            }
        }
#endif /* FOAM_ADIOS_CLOUD_EXPAND */

        outputSize_ += myParticles *
        (
            nLabels    * adiosTraits<label>::adiosSize
          + nScalars   * adiosTraits<scalar>::adiosSize
          + 3*nVectors * adiosTraits<scalar>::adiosSize
        );
    }

    if (nClouds)
    {
        const fileName varName = r.regionPath();

        // number of active clouds as region attribute
        defineIntAttribute("nClouds", varName, nClouds);
    }

    return maxLen;
}


void Foam::adiosWrite::cloudWrite(const regionInfo& r)
{
    const fvMesh& mesh = time_.lookupObject<fvMesh>(r.name_);
    Info<< "  adiosWrite::cloudWrite: region "
        << r.index_ << " " << r.name_ << ": " << endl;

#ifdef FOAM_ADIOS_CLOUD_EXPLICIT_NAMES
    cloudWriteExplicit(r);
#endif

    DynamicList<label> labelBuffer;
    DynamicList<scalar> scalarBuffer;

    label nClouds = 0;
    forAll(r.cloudNames_, cloudI)
    {
        Info<< "    cloud: " << r.cloudNames_[cloudI] << endl;
        Info<< "    Properties " <<  basicKinematicCloud::particleType::propertyList() << endl;

        // If the cloud contains no particles, jump to the next cloud
        if (r.nTotalParticles_ == 0)
        {
            continue;
        }

        const kinematicCloud& cloud =
            mesh.lookupObject<kinematicCloud>(r.cloudNames_[cloudI]);
        const basicKinematicCloud& q =
            static_cast<const basicKinematicCloud&>(cloud);

        // Number of particles on this process
        const label myParticles = q.size();

        const fileName varPath = r.cloudPath(nClouds++);

        // number of particles per processor as a field
        writeIntVariable(varPath/"nParticle", myParticles);

        // stream contents
        {
            ORawBufStream os(iobuffer_);  // always raw binary content

            forAllConstIter(basicKinematicCloud, q, pIter)
            {
                os << *pIter;
            }
        }

        writeVariable(varPath/"__blob__", iobuffer_);

        labelBuffer.reserve(myParticles);
        scalarBuffer.reserve(3*myParticles);

        // buffer for 'label'
        labelBuffer.setSize(myParticles);

        // buffer for 'vector'
        scalarBuffer.setSize(3*myParticles);

        // Write slip velocity Us = U - Uc
        {
            const word what = "Us";
            if (findStrings(r.cloudAttribs_, what))
            {
                label i = 0;
                forAllConstIter(basicKinematicCloud, q, pIter)
                {
                    scalarBuffer[3*i+0] = pIter().U().x() - pIter().Uc().x();
                    scalarBuffer[3*i+1] = pIter().U().y() - pIter().Uc().y();
                    scalarBuffer[3*i+2] = pIter().U().z() - pIter().Uc().z();
                    ++i;
                }
                writeVariable(varPath/what, scalarBuffer);
            }
        }

#ifdef FOAM_ADIOS_CLOUD_EXPAND

        // expanding particle blob into separate fields
        // mostly useful for debugging and as a general example of working
        // with particle blobs

        // walk the blob framents
        labelBuffer.reserve(myParticles);
        scalarBuffer.reserve(3*myParticles);

        ParticleBinaryBlob blob
        (
            basicKinematicCloud::particleType::propertyTypes(),
            basicKinematicCloud::particleType::propertyList()
        );

        List<char> blobBuffer(blob.byteSize() + 32); // extra generous
        ORawBufStream osblob(blobBuffer);

        labelBuffer.setSize(myParticles);
        scalarBuffer.setSize(3*myParticles);

        forAllConstIter(SLList<ParticleBinaryBlob::Fragment>, blob, iter)
        {
            const ParticleBinaryBlob::Fragment& frag = *iter;

            Info<< frag << endl;

            if (frag.type() == pTraits<label>::typeName)
            {
                label i = 0;
                forAllConstIter(basicKinematicCloud, q, pIter)
                {
                    osblob.rewind();
                    osblob << *pIter;   // binary content

                    labelBuffer[i] = frag.getLabel(blobBuffer);
                    ++i;
                }

                writeVariable(varPath/frag.name(), labelBuffer);
            }
            else if (frag.type() == pTraits<scalar>::typeName)
            {
                label i = 0;
                forAllConstIter(basicKinematicCloud, q, pIter)
                {
                    osblob.rewind();
                    osblob << *pIter;   // binary content

                    scalarBuffer[i] = frag.getScalar(blobBuffer);
                    ++i;
                }

                writeVariable(varPath/frag.name(), scalarBuffer);
            }
            else if (frag.type() == pTraits<vector>::typeName)
            {
                label i = 0;
                forAllConstIter(basicKinematicCloud, q, pIter)
                {
                    osblob.rewind();
                    osblob << *pIter;   // binary content
                    vector val = frag.getVector(blobBuffer);

                    scalarBuffer[3*i+0] = val.x();
                    scalarBuffer[3*i+1] = val.y();
                    scalarBuffer[3*i+2] = val.z();
                    ++i;
                }

                writeVariable(varPath/frag.name(), scalarBuffer);
            }
        }

#endif /* FOAM_ADIOS_CLOUD_EXPAND */
    }

}


size_t Foam::adiosWrite::cloudDefineExplicit(regionInfo& r)
{
#ifdef FOAM_ADIOS_CLOUD_EXPLICIT_NAMES
    Info<< "  adiosWrite::cloudDefineExplicit: region "
        << r.index_ << " " << r.name_ << ": " << endl;

    const fvMesh& mesh = time_.lookupObject<fvMesh>(r.name_);

    forAll(r.cloudNames_, cloudI)
    {
        Info<< "    cloud: " << r.cloudNames_[cloudI] << endl;

        const kinematicCloud& cloud =
            mesh.lookupObject<kinematicCloud>(r.cloudNames_[cloudI]);

        const basicKinematicCloud& q =
            static_cast<const basicKinematicCloud&>(cloud);

        // Number of particles on this process
        label myParticles = q.size();

        // Find the number of particles on each process
        r.nParticles_[Pstream::myProcNo()] = myParticles;
        Pstream::gatherList(r.nParticles_);
        Pstream::scatterList(r.nParticles_);

        // Sum total number of particles on all processes
        r.nTotalParticles_ = sum(r.nParticles_);

        // If the cloud contains no particles, jump to the next cloud
        if (r.nTotalParticles_ == 0)
        {
            Info<< "    " << r.cloudNames_[cloudI]
                <<": No particles in cloud. Skipping definition."
                << endl;

            continue;
        }

        // Calculate offset values
        List<label> offsets(Pstream::nProcs());
        offsets[0] = 0;
        for (label proc=1; proc < offsets.size(); proc++)
        {
            offsets[proc] = offsets[proc-1] + r.nParticles_[proc-1];
        }

        // Define all possible output variables, not necessary to write all of them later
        const fileName varPath = oldCloudPath(r.index_, r.cloudNames_[cloudI]);

        defineIntVariable(varPath/"nParticlesPerProc");
        defineIntVariable(varPath/"nTotalParticles");

        // Define a variable for dataset name
        char gdimstr[16]; // 1D global array of particles from all processes
        char ldimstr[16]; // this process' size
        char offsstr[16]; // this process' offset in global array

        sprintf(gdimstr, "%d", r.nTotalParticles_);
        sprintf(ldimstr, "%d", myParticles);
        sprintf(offsstr, "%d", offsets[Pstream::myProcNo()]);

        // integer info datasets
        adios_define_var(groupID_, (varPath/"origProc").c_str(), "", adios_integer, ldimstr, gdimstr, offsstr);
        adios_define_var(groupID_, (varPath/"origId").c_str(),   "", adios_integer, ldimstr, gdimstr, offsstr);
        adios_define_var(groupID_, (varPath/"cell").c_str(),     "", adios_integer, ldimstr, gdimstr, offsstr);
        adios_define_var(groupID_, (varPath/"currProc").c_str(), "", adios_integer, ldimstr, gdimstr, offsstr);
        label nLabels = 4;

        // scalar datasets
        adios_define_var(groupID_, (varPath/"rho").c_str(), "", ADIOS_IOSCALAR, ldimstr, gdimstr, offsstr);
        adios_define_var(groupID_, (varPath/"d").c_str(),   "", ADIOS_IOSCALAR, ldimstr, gdimstr, offsstr);
        adios_define_var(groupID_, (varPath/"age").c_str(), "", ADIOS_IOSCALAR, ldimstr, gdimstr, offsstr);
        label nScalars = 3;

        // vector datasets
        sprintf(gdimstr, "%d,3", r.nTotalParticles_);
        sprintf(ldimstr, "%d,3", myParticles);
        sprintf(offsstr, "%d",   offsets[Pstream::myProcNo()]);

        adios_define_var(groupID_, (varPath/"position").c_str(), "", ADIOS_IOSCALAR, ldimstr, gdimstr, offsstr);
        adios_define_var(groupID_, (varPath/"U").c_str(),        "", ADIOS_IOSCALAR, ldimstr, gdimstr, offsstr);
        adios_define_var(groupID_, (varPath/"Us").c_str(),       "", ADIOS_IOSCALAR, ldimstr, gdimstr, offsstr);
        label nVectors = 3;

        outputSize_ += myParticles *
        (
            nLabels    * adiosTraits<label>::adiosSize
          + nScalars   * sizeof(ioScalar)
          + 3*nVectors * sizeof(ioScalar)
        );
    }
#endif /* FOAM_ADIOS_CLOUD_EXPLICIT_NAMES */

    return 0;
}


void Foam::adiosWrite::cloudWriteExplicit(const regionInfo& r)
{
#ifdef FOAM_ADIOS_CLOUD_EXPLICIT_NAMES

    DynamicList<label> labelBuffer;
    DynamicList<scalar> scalarBuffer;

    const fvMesh& mesh = time_.lookupObject<fvMesh>(r.name_);
    Info<< "  adiosWrite::cloudWriteExplicit: region "
        << r.index_ << " " << r.name_ << ": " << endl;

    forAll(r.cloudNames_, cloudI)
    {
        Info<< "    cloud: " << r.cloudNames_[cloudI] << endl;
        Info<< "    Properties " <<  basicKinematicCloud::particleType::propertyList() << endl;

        // If the cloud contains no particles, jump to the next cloud
        if (r.nTotalParticles_ == 0)
        {
            continue;
        }

        const kinematicCloud& cloud =
            mesh.lookupObject<kinematicCloud>(r.cloudNames_[cloudI]);
        const basicKinematicCloud& q =
            static_cast<const basicKinematicCloud&>(cloud);

        const label myParticles = q.size();
        const fileName varPath = oldCloudPath(r.index_, r.cloudNames_[cloudI]);

        writeIntVariable(varPath/"nParticlesPerProc", myParticles);
        writeIntVariable(varPath/"nTotalParticles",   r.nTotalParticles_);

        labelBuffer.reserve(myParticles);
        scalarBuffer.reserve(3*myParticles);

        // buffer for 'label' and scalar/vector
        labelBuffer.setSize(myParticles);
        scalarBuffer.setSize(3*myParticles);

        // Write original processor ID
        {
            const word what = "origProc";
            if (findStrings(r.cloudAttribs_, what))
            {
                Info<< "      dataset " << what << endl;
                label i = 0;
                forAllConstIter(basicKinematicCloud, q, pIter)
                {
                    labelBuffer[i++] = pIter().origProc();
                }
                writeVariable(varPath/what, labelBuffer);
            }
        }

        // Write original ID
        {
            const word what = "origId";
            if (findStrings(r.cloudAttribs_, what))
            {
                Info<< "      dataset " << what << endl;
                label i = 0;
                forAllConstIter(basicKinematicCloud, q, pIter)
                {
                    labelBuffer[i++] = pIter().origId();
                }
                writeVariable(varPath/what, labelBuffer);
            }
        }

        // Write cell number
        {
            const word what = "cell";
            if (findStrings(r.cloudAttribs_, what))
            {
                Info<< "      dataset " << what << endl;
                label i = 0;
                forAllConstIter(basicKinematicCloud, q, pIter)
                {
                    labelBuffer[i++] = pIter().cell();
                }
                writeVariable(varPath/what, labelBuffer);
            }
        }

        // Write current process ID
        {
            const word what = "currProc";
            if (findStrings(r.cloudAttribs_, what))
            {
                Info<< "      dataset " << what << endl;
                label i = 0;
                forAllConstIter(basicKinematicCloud, q, pIter)
                {
                    labelBuffer[i++] = Pstream::myProcNo();
                }
                writeVariable(varPath/what, labelBuffer);
            }
        }

        // Write density rho
        {
            const word what = "rho";
            if (findStrings(r.cloudAttribs_, "rho"))
            {
                Info<< "      dataset " << what << endl;
                label i = 0;
                forAllConstIter(basicKinematicCloud, q, pIter)
                {
                    scalarBuffer[i++] = pIter().rho();
                }
                writeVariable(varPath/what, scalarBuffer);
            }
        }

        // Write diameter d
        {
            const word what = "d";
            if (findStrings(r.cloudAttribs_, what))
            {
                Info<< "      dataset " << what << endl;
                label i = 0;
                forAllConstIter(basicKinematicCloud, q, pIter)
                {
                    scalarBuffer[i++] = pIter().d();
                }
                writeVariable(varPath/what, scalarBuffer);
            }
        }

        // Write age
        {
            const word what = "age";
            if (findStrings(r.cloudAttribs_, what))
            {
                Info<< "      dataset " << what << endl;
                label i = 0;
                forAllConstIter(basicKinematicCloud, q, pIter)
                {
                    scalarBuffer[i++] = pIter().age();
                }
                writeVariable(varPath/what, scalarBuffer);
            }
        }

        // buffer for 'vector'
        scalarBuffer.setSize(3*myParticles);

        // Write position
        {
            const word what = "position";
            if (findStrings(r.cloudAttribs_, what))
            {
                label i = 0;
                forAllConstIter(basicKinematicCloud, q, pIter)
                {
                    scalarBuffer[3*i+0] = pIter().position().x();
                    scalarBuffer[3*i+1] = pIter().position().y();
                    scalarBuffer[3*i+2] = pIter().position().z();
                    ++i;
                }
                writeVariable(varPath/what, scalarBuffer);
            }
        }

        // Write velocity U
        {
            const word what = "U";
            if (findStrings(r.cloudAttribs_, what))
            {
                label i = 0;
                forAllConstIter(basicKinematicCloud, q, pIter)
                {
                    scalarBuffer[3*i+0] = pIter().U().x();
                    scalarBuffer[3*i+1] = pIter().U().y();
                    scalarBuffer[3*i+2] = pIter().U().z();
                    ++i;
                }
                writeVariable(varPath/what, scalarBuffer);
            }
        }

        // Write slip velocity Us = U - Uc
        {
            const word what = "Us";
            if (findStrings(r.cloudAttribs_, what))
            {
                label i = 0;
                forAllConstIter(basicKinematicCloud, q, pIter)
                {
                    scalarBuffer[3*i+0] = pIter().U().x() - pIter().Uc().x();
                    scalarBuffer[3*i+1] = pIter().U().y() - pIter().Uc().y();
                    scalarBuffer[3*i+2] = pIter().U().z() - pIter().Uc().z();
                    ++i;
                }
                writeVariable(varPath/what, scalarBuffer);
            }
        }
    }

#endif /* FOAM_ADIOS_CLOUD_EXPLICIT_NAMES */
}


// ************************************************************************* //
