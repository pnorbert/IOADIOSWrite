/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright(C) 2015 Norbert Podhorszki
     \\/     M anipulation  | Copyright(C) 2016 OpenCFD Ltd.
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
#include "basicKinematicCollidingCloud.H"
#include "reactingCloud.H"
#include "emptyFvPatchField.H"

#include "IBufStream.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::adiosTime Foam::adiosWrite::readData(const fileName& bpFile)
{
    Info<< " Read data of step " << bpFile << endl;

    adiosReader reader(bpFile.c_str(), comm_);

    if (!reader.isGood())
    {
        return adiosTime();
    }

    adiosTime timeInfo(reader);

    iobuffer_.reserve(reader.sizeOf());

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

        HashTable<adiosReader::fieldInfo> fromFile = reader.getFieldInfo(regName);
        wordList fieldNames = fromFile.sortedToc();

        forAll(fieldNames, fieldI)
        {
            const word& name = fieldNames[fieldI];
            const adiosReader::fieldInfo& src = fromFile[name];

            if ((static_cast<objectRegistry&>(mesh)).found(name))
            {
                readVolField(mesh.find(name)(), reader, src);
            }
        }
    }


    // clouds
    forAll(meshNames, regI)
    {
        const word& regName = meshNames[regI];
        fvMesh& mesh = const_cast<fvMesh&>(time_.lookupObject<fvMesh>(regName));

        Info<<"lookup: " << regName << endl;

        HashTable<adiosReader::cloudInfo> fromFile = reader.getCloudInfo(regName);
        wordList cloudNames = fromFile.sortedToc();

        forAll(cloudNames, fieldI)
        {
            const word& name = cloudNames[fieldI];
            const adiosReader::cloudInfo& src = fromFile[name];

            if ((static_cast<objectRegistry&>(mesh)).found(name))
            {
                Info<<"read-back cloud " << name << endl;
                readCloud(mesh, name, reader, src);
            }
         }
    }

    reader.close();

    return timeInfo;
}


Foam::adiosTime Foam::adiosWrite::readData(const instant& when)
{
    Info<< " Read data of step " << when.name() << endl;
    return readData(when.name());
}


Foam::adiosTime Foam::adiosWrite::readData()
{
    Info<< " Read data of step " << obr_.time().timeName() << endl;
    return readData(dataDirectory/obr_.time().timeName() + ".bp");
}


bool Foam::adiosWrite::readVolField
(
    regIOobject* obj,
    const adiosReader& reader,
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
            reader,
            src
        );
    }
    else if (fieldType == volVectorField::typeName)
    {
        return fieldRead
        (
            static_cast<volVectorField&>(*obj),
            reader,
            src
        );
    }
    else if (fieldType == surfaceScalarField::typeName)
    {
        return fieldRead
        (
            static_cast<surfaceScalarField&>(*obj),
            reader,
            src
        );
    }
    else if (fieldType == volSphericalTensorField::typeName)
    {
        return fieldRead
        (
            static_cast<volSphericalTensorField&>(*obj),
            reader,
            src
        );
    }
    else if (fieldType == volSymmTensorField::typeName)
    {
        return fieldRead
        (
            static_cast<volSymmTensorField&>(*obj),
            reader,
            src
        );
    }
    else if (fieldType == volTensorField::typeName)
    {
        return fieldRead
        (
            static_cast<volTensorField&>(*obj),
            reader,
            src
        );
    }
    else
    {
        return false;
    }
}


bool Foam::adiosWrite::readCloud
(
    const fvMesh& mesh,
    const word& cloudName,
    const adiosReader& reader,
    const adiosReader::cloudInfo& src
)
{
    regIOobject* obj = *(const_cast<fvMesh&>(mesh).find(cloudName));
    const word& cloudType = obj->type();
    const word&   srcType = src.type();

    if (cloudType != srcType)
    {
        // probably fatal:
        Info<<"WARNING mismatch on cloud " << src
            << " (expected: " << cloudType
            << " but had " << srcType << ")\n";
        return false;
    }

    if (cloudType == Cloud<basicKinematicCollidingParcel>::typeName)
    {
        Info<<"read cloud " << obj->name() << " (type " << cloudType
            << ") from file\n";

#if 0
        // this is just to check the round-trip information
        // retrieve particle fragments
        if (src.nParticle())
        {
            const size_t incr = src.width();
            const ParticleBinaryBlob& binfo = src.blob();

            size_t nread = reader.getBuffered(src.fullName(), iobuffer_);
            char* rawdata = reinterpret_cast<char*>(iobuffer_.data());

            for (int pI=0; pI < src.nParticle(); ++pI)
            {
                const char *blobBuffer = (rawdata + pI*incr);

                forAllConstIter(ParticleBinaryBlob::Container, binfo, iter)
                {
                    const ParticleBinaryBlob::Fragment& frag = *iter;

                    if (frag.type() == pTraits<label>::typeName)
                    {
                        Info<< frag.name() << " (" << frag.type() << ") = "
                            << frag.getLabel(blobBuffer)  << nl;
                    }
                    else if (frag.type() == pTraits<scalar>::typeName)
                    {
                        Info<< frag.name() << " (" << frag.type() << ") = "
                            << frag.getScalar(blobBuffer)  << nl;
                    }
                    else if (frag.type() == pTraits<vector>::typeName)
                    {
                        Info<< frag.name() << " (" << frag.type() << ") = "
                            << frag.getVector(blobBuffer)  << nl;
                    }
                }
            }

            // not really serious about doing this, but for diagnostics

            OStringStream os(IOstream::ASCII);

            os << src.nParticle() << '(';

            for (int pI=0; pI < src.nParticle(); ++pI)
            {
                const char *blobBuffer = (rawdata + pI*incr);

                forAllConstIter(ParticleBinaryBlob::Container, binfo, iter)
                {
                    const ParticleBinaryBlob::Fragment& frag = *iter;

                    if (frag.type() == pTraits<label>::typeName)
                    {
                        os << ' ' << frag.getLabel(blobBuffer);
                    }
                    else if (frag.type() == pTraits<scalar>::typeName)
                    {
                        os << ' ' << frag.getScalar(blobBuffer);
                    }
                    else if (frag.type() == pTraits<vector>::typeName)
                    {
                        os << ' ' << frag.getVector(blobBuffer);
                    }
                }
            }

            os << ')';

        }
#endif

        if (true)
        {
            Info<<"skipping .. currently not working properly\n";
            return false;
        }

        // TODO: other cloud types
        typedef basicKinematicCloud::particleType pType;
        basicKinematicCloud& cloudObj = static_cast<basicKinematicCloud&>(*obj);

        Pout<<"Current cloud: " << cloudObj.size() << " parcels. Cloud on file: "
            << src.nParticle() << " elements" << endl;

        IDLList<pType> newParticles;

        if (cloudObj.size() > src.nParticle())
        {
            // TODO: truncate!
            Pout<<"WARNING cloud " << obj->name() << " needs truncating. Has "
                << cloudObj.size() << " parcels, but file has "
                << src.nParticle() << endl;
        }
        else if (cloudObj.size() < src.nParticle())
        {
            // TODO: add to end of list!
            Pout<<"WARNING cloud " << obj->name() << " needs expansion. Has "
                << cloudObj.size() << " parcels, but file has "
                << src.nParticle() << endl;
        }

        if (src.nParticle())
        {
            const size_t incr = src.width();
            size_t nread = reader.getBuffered(src.fullName(), iobuffer_);

#if 1
            const ParticleBinaryBlob& binfo = src.blob();
            char* rawdata = reinterpret_cast<char*>(iobuffer_.data());

            // not really serious about doing this, but perhaps usable for debugging

            OStringStream os(IOstream::ASCII);
            os << src.nParticle() << '(';

            for (int pI=0; pI < src.nParticle(); ++pI)
            {
                const char *blobBuffer = (rawdata + pI*incr);

                forAllConstIter(ParticleBinaryBlob::Container, binfo, iter)
                {
                    const ParticleBinaryBlob::Fragment& frag = *iter;

                    if (frag.type() == pTraits<label>::typeName)
                    {
                        os << ' ' << frag.getLabel(blobBuffer);
                    }
                    else if (frag.type() == pTraits<scalar>::typeName)
                    {
                        os << ' ' << frag.getScalar(blobBuffer);
                    }
                    else if (frag.type() == pTraits<vector>::typeName)
                    {
                        os << ' ' << frag.getVector(blobBuffer);
                    }
                    else
                    {
                        // Fatal: frag.type()
                    }
                }
            }
            os << ')';

            IStringStream is(os.str(), IOstream::ASCII);

            newParticles = IDLList<pType>
            (
                is,
                pType::iNew(mesh)
            );

#else

            // one particle
            label num = 1;

            char databuf[num*(incr + 4) + 4];

            char *src = reinterpret_cast<char*>(iobuffer_.data());
            char *dst = databuf;
            size_t nbytes = 0;

            dst[nbytes++] = '(';
            for (int pI=0; pI < num; ++pI)
            {
                dst[nbytes++] = '(';

                for (size_t i=0; i < incr; ++i)
                {
                    dst[nbytes++] = *(src++);
                }

                dst[nbytes++] = ')';
            }
            dst[nbytes++] = ')';


            string input(databuf, nbytes);
            Pout<< "input (" << nbytes << ") " << input << endl;

            // testing - basic transcription:
            IBufStream is(databuf, nbytes, IOstream::BINARY);
            newParticles = IDLList<pType>
            (
                is,
                pType::iNew(mesh)
            );
#endif
///            {
///                nbytes += src.sizeOf();
///            }

///            label ndigits = 1 + floor(Foam::log10(scalar(src.nParticle())));
///            // include surrounding '()'
///            size_t nbytes = src.sizeOf() + ndigits + 2;
///
///            const size_t sizeT = sizeof(DynamicCharList::value_type);
///            const size_t target = (nbytes / sizeT) + (nbytes % sizeT);
///
///            iobuffer_.setSize(target);
///            char *data = reinterpret_cast<char*>(iobuffer_.data());
///
///            OBufStream os(data, nbytes, IOstream::ASCII);
///
///            os << src.nParticle() << '(';
///            nbytes = os.size();
///
///            if (reader.getVariable(src.fullName(), (data+nbytes)))
///            {
///                nbytes += src.sizeOf();
///            }
///            data[nbytes] = ')';
///
///            string input(data, nbytes);
///            Pout<< "Got input (" << nbytes << ") " << input << endl;
///
///            IBufStream is(data, nbytes, IOstream::BINARY);
///            newParticles = IDLList<pType>
///            (
///                is,
///                pType::iNew(mesh)
///            );
        }

#if 1
        cloudObj.clear();

        forAllIter(Cloud<pType>, newParticles, newpIter)
        {
            pType& newp = newpIter();

            cloudObj.addParticle(newParticles.remove(&newp));
        }
#endif
        return true;
    }
    else
    {
        Info<<"read cloud " << obj->name() << " (type " << cloudType
            << ") from file - unsupported type" << endl;


        return false;
    }
}


// ************************************************************************* //
