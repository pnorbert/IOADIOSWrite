/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
     \\/     M anipulation  |
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

template<class CloudType>
size_t Foam::adiosWrite::cloudDefine
(
    const CloudType& cld,
    cloudInfo& cInfo,
    const fileName& varPath
)
{
    typedef OCompactCountStream ParticleCountStream;

    size_t bufLen = 0;
    size_t maxLen = 0;

    const word& cloudName = cInfo.name();
    const word& cloudType = cInfo.type();

    // Number of particles on this process
    const label nParticle = cld.size();

    // Set number of particles on process and total of all processes
    if (cInfo.nParticle(nParticle) == 0)
    {
        // skip: cloud has no particles
        Info<< "    " << cloudName
            <<": No particles in cloud. Skipping definition."
            << endl;

        return 0;
    }

    //
    // determine binary representation
    //

    // from compile-time information
    ParticleBinaryBlob binfo
    (
        CloudType::particleType::propertyTypes(),
        CloudType::particleType::propertyList()
    );


    // from run-time information
    List<label> blobSize(Pstream::nProcs(), 0);
    if (nParticle)
    {
        ParticleCountStream os(IOstream::BINARY);  // binary content
        os << *(cld.first());

        blobSize[Pstream::myProcNo()] = os.size();
    }
    Pstream::gatherList(blobSize);
    Pstream::scatterList(blobSize);

    if (Foam::max(blobSize) != label(binfo.byteSize()))
    {
        // probably also fatal:
        WarningInFunction
            << "Mismatch in blob-sizes for cloud " << cloudName
            << " expected blob-size: " <<  binfo.byteSize()
            << " found blob-size: " << blobSize << nl
            << "  types: " << FlatListOutput<word>(binfo.types()) << nl
            << "  names: " << FlatListOutput<word>(binfo.names()) << nl
            << "  offset:" << FlatListOutput<label>(binfo.offsets()) << nl
            << "  width: " << FlatListOutput<label>(binfo.sizes()) << endl;

        return 0;
    }

    // added verbosity
    if (1)
    {
        Info<< "cloud: " << cloudName << " (" << cInfo.nTotal()
            << " particles, " << binfo.byteSize() << " bytes) type: "
            << cloudType << nl
            << "  types: " << FlatListOutput<word>(binfo.types()) << nl
            << "  names: " << FlatListOutput<word>(binfo.names()) << nl
            << "  offset:" << FlatListOutput<label>(binfo.offsets()) << nl
            << "  width: " << FlatListOutput<label>(binfo.sizes()) << endl;
    }
    else
    {
        Info<< "cloud: " << cloudName << " (" << cInfo.nTotal()
            << " particles, " << binfo.byteSize() << " bytes) type: "
            << cloudType << nl;
    }


    //
    // cloud attributes:
    //

    defineAttribute("class", varPath, cloudType);

    // total number of particles as an attribute
    defineIntAttribute("nParticle",  varPath, cInfo.nTotal());

    defineIntAttribute("size",       varPath, binfo.byteSize());
    defineListAttribute("names",     varPath, binfo.names());
    defineListAttribute("types",     varPath, binfo.types());
    defineListAttribute("offset",    varPath, binfo.offsets());
    defineListAttribute("byte-size", varPath, binfo.sizes());


    //
    // cloud variables:
    //

    // number of particles per processor as a variable
    defineIntVariable(varPath/"nParticle");


    // common - local offset value
    const string offsetDims = Foam::name(cInfo.offset());

    // for scalars
    const string localDims  = Foam::name(nParticle);
    const string globalDims = Foam::name(cInfo.nTotal());

    // for vectors
    const string localDims3  = localDims  + ",3";
    const string globalDims3 = globalDims + ",3";

    {
        // stream contents
        const string localBlob  = localDims  + "," + Foam::name(binfo.byteSize());
        const string globalBlob = globalDims + "," + Foam::name(binfo.byteSize());

        adios_define_var
        (
            groupID_,
            varPath.c_str(),                // name
            NULL,                           // path (deprecated)
            adios_unsigned_byte,            // data-type
            localBlob.c_str(),              // local dimensions
            globalBlob.c_str(),             // global dimensions
            offsetDims.c_str()              // local offsets
        );

        bufLen = nParticle * binfo.byteSize();
        maxLen = Foam::max(maxLen, bufLen);

        outputSize_ += bufLen;
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
        forAllConstIter(ParticleBinaryBlob::Container, binfo, iter)
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

    outputSize_ += nParticle *
    (
        nLabels    * adiosTraits<label>::adiosSize
      + nScalars   * adiosTraits<scalar>::adiosSize
      + 3*nVectors * adiosTraits<scalar>::adiosSize
    );

    return maxLen;
}


template<class CloudType>
void Foam::adiosWrite::cloudWrite
(
    const CloudType& cld,
    const cloudInfo& cInfo,
    const regionInfo& rInfo
)
{
    typedef OCompactBufStream ParticleOutputStream;

    DynamicList<label> labelBuffer;
    DynamicList<scalar> scalarBuffer;

    const word& cloudName = cInfo.name();
    const fileName varPath = rInfo.cloudPath(cloudName);


    // If the cloud contains no particles, jump to the next cloud
    // - this is probably redundant here
    if (cInfo.nTotal() == 0)
    {
        return;
    }

    Info<< "    cloud: " << cloudName << endl;

    // Number of particles on this process
    const label nParticle = cld.size();

    // number of particles per processor as a field
    writeIntVariable(varPath/"nParticle", nParticle);

    // TODO: reserve enough space on iobuffer

    // stream the cloud contents - always binary content
    {
        ParticleOutputStream os(iobuffer_, IOstream::BINARY);

        forAllConstIter(typename CloudType, cld, pIter)
        {
            os << *pIter;
        }
    }

    writeVariable(varPath, iobuffer_);

    labelBuffer.reserve(nParticle);
    scalarBuffer.reserve(3*nParticle);

    // buffer for 'label'
    labelBuffer.setSize(nParticle);

    // buffer for 'vector'
    scalarBuffer.setSize(3*nParticle);

    // Write slip velocity Us = U - Uc

    // Cannot be in the template!
    // eg, error: 'const class Foam::indexedParticle' has no member named 'Uc'

#if 0
    {
        const word what = "Us";
        if (findStrings(rInfo.requestedAttrs_, what))
        {
            label i = 0;
            forAllConstIter(typename CloudType, cld, pIter)
            {
                scalarBuffer[3*i+0] = pIter().U().x() - pIter().Uc().x();
                scalarBuffer[3*i+1] = pIter().U().y() - pIter().Uc().y();
                scalarBuffer[3*i+2] = pIter().U().z() - pIter().Uc().z();
                ++i;
            }
            writeVariable(varPath/what, scalarBuffer);
        }
    }
#endif

#ifdef FOAM_ADIOS_CLOUD_EXPAND

    // expanding particle blob into separate fields
    // mostly useful for debugging and as a general example of working
    // with particle blobs

    // walk the blob framents
    labelBuffer.reserve(nParticle);
    scalarBuffer.reserve(3*nParticle);

    ParticleBinaryBlob binfo
    (
        CloudType::particleType::propertyTypes(),
        CloudType::particleType::propertyList()
    );

    List<char> blobBuffer(binfo.byteSize() + 32); // extra generous
    ParticleOutputStream osblob(blobBuffer, IOstream::BINARY);

    labelBuffer.setSize(nParticle);
    scalarBuffer.setSize(3*nParticle);

    forAllConstIter(ParticleBinaryBlob::Container, binfo, iter)
    {
        const ParticleBinaryBlob::Fragment& frag = *iter;

        if (frag.type() == pTraits<label>::typeName)
        {
            label i = 0;
            forAllConstIter(typename CloudType, cld, pIter)
            {
                osblob.rewind();
                osblob << *pIter;   // binary content

                labelBuffer[i] = frag.getLabel(blobBuffer);
                ++i;
            }

            // Info<< frag.name() << " (" << frag.type() << ") = "
            //     << FlatListOutput<int>(SubList<int>(labelBuffer, i)) << nl;

            writeVariable(varPath/frag.name(), labelBuffer);
        }
        else if (frag.type() == pTraits<scalar>::typeName)
        {
            label i = 0;
            forAllConstIter(typename CloudType, cld, pIter)
            {
                osblob.rewind();
                osblob << *pIter;   // binary content

                scalarBuffer[i] = frag.getScalar(blobBuffer);
                ++i;
            }

            // Info<< frag.name() << " (" << frag.type() << ") = "
            //     << FlatListOutput<scalar>(SubList<scalar>(scalarBuffer, i)) << nl;

            writeVariable(varPath/frag.name(), scalarBuffer);
        }
        else if (frag.type() == pTraits<vector>::typeName)
        {
            label i = 0;
            forAllConstIter(typename CloudType, cld, pIter)
            {
                osblob.rewind();
                osblob << *pIter;   // binary content
                vector val = frag.getVector(blobBuffer);

                scalarBuffer[3*i+0] = val.x();
                scalarBuffer[3*i+1] = val.y();
                scalarBuffer[3*i+2] = val.z();
                ++i;
            }

            // Info<< frag.name() << " (" << frag.type() << ") = "
            //     << FlatListOutput<scalar>(SubList<scalar>(scalarBuffer, i)) << nl;

            writeVariable(varPath/frag.name(), scalarBuffer);
        }
    }

#endif /* FOAM_ADIOS_CLOUD_EXPAND */
}


// ************************************************************************* //
