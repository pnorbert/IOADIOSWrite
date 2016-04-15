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
#include "IStringStream.H"
#include "IStringStreamBuf.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::adiosWrite::readData(const fileName& bpFile)
{
    Info<< " Read data of step " << bpFile << endl;

    adiosReader::helper helper;

    if (!helper.open(bpFile.c_str(), comm_))
    {
        return false;
    }

    // Process N will read writeblock N
    helper.select(adios_selection_writeblock(Pstream::myProcNo()));

    helper.scan(true);
    IStringStreamBuf is(helper.maxLen, IOstream::BINARY);

    bool ok = true;
    forAll(regions_, regionId)
    {
        regionInfo& rInfo = regions_[regionId];

        ok =
        (
            ok
         && fieldRead(is, helper, rInfo)
         && readClouds(helper, rInfo)
        );
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


bool Foam::adiosWrite::fieldRead
(
    IStringStreamBuf& is,
    adiosReader::helper& helper,
    regionInfo& rInfo
)
{
    const fvMesh& mesh = time_.lookupObject<fvMesh>(rInfo.name_);

    Pout<< "  adiosWrite::fieldRead: region "
        << rInfo.index_ << "=" << rInfo.name_ << endl;

    return
    (
        fieldRead<volScalarField>
        (
            is,
            helper,
            mesh,
            rInfo,
            rInfo.scalarFields_
        )
     && fieldRead<volVectorField>
        (
            is,
            helper,
            mesh,
            rInfo,
            rInfo.vectorFields_
        )
     && fieldRead<surfaceScalarField>
        (
            is,
            helper,
            mesh,
            rInfo,
            rInfo.surfaceScalarFields_
        )
    );
}



// processor0/0.004/uniform/lagrangian/kinematicCloud/cloudProperties <stream>   - same for all processors
// processor0/0.004/uniform/lagrangian/kinematicCloud/kinematicCloudOutputProperties <stream>    - same for all processors

// processor0/0.004/lagrangian/kinematicCloud/
// processor0/0.004/lagrangian/kinematicCloud/{U, age, active, angularMom} ...

/// bool Foam::adiosWrite::readScalarFields
/// (
///     IStringStreamBuf& is,
///     adiosReader::helper& helper,
///     regionInfo& rInfo
/// )
/// {
///     bool ok = true;
///
///     const fvMesh& mesh = time_.lookupObject<fvMesh>(regions_[regionID].name_);
/// #if 1
///     ok = fieldRead<volScalarField>
///     (
///         is,
///         helper,
///         mesh,
///         regions_[regionID].scalarFields_,
///         regionID
///     );
/// #else
///     const fieldGroup<scalar>& fields(regions_[regionID].scalarFields_);
///
///     forAll(fields, fieldI)
///     {
///         // Lookup field
///         volScalarField& field =
///             const_cast<volScalarField&>
///             (
///                 mesh.lookupObject<volScalarField>(fields[fieldI])
///             );
///
///         Info<< "    readScalarField via dictionary: " << field.name() << endl;
///
///         fileName datasetName = rInfo.fieldVarPath(fields[fieldI]);
///
///         // is.reserve(sizeOf(datasetName, true));
///         // return getDataSet(datasetName, is.data());
///
///         // Read data from file
///         ok = helper.getDataSet(datasetName, is);
///         if (ok)
///         {
///             Info<<"istream content:" << endl;
///             Info<<"has " << is.stdStream().rdbuf()->in_avail() << " chars" << endl;
///             is.print(Info);
///             Info<< "is.good: " << is.good() << endl;
///             Info<< "is.eof: " << is.eof() << endl;
///
///             char c;
///             while (is.good() && !is.eof())
///             {
///                 is.get(c);
///                 Info<< char(c);
///             }
///             Info<< endl << "DONE" << endl;
///
///             // read fields via dictionary
///             // dictionary dict(is);
///
///             // Info<<"dictionary: " << dict << endl;
///
///             // field.readField(dict, "internalField");
///             // field.boundaryField().readField(field, dict.subDict("boundaryField"));
///
///             // TODO: adjust for referenceLevel?
///             //
///             //     if (dict.found("referenceLevel"))
///             //     {
///             //         Type fieldAverage(pTraits<Type>(dict.lookup("referenceLevel")));
///             //         Field<Type>::operator+=(fieldAverage);
///             //         forAll(boundaryField_, patchI)
///             //         {
///             //             boundaryField_[patchI] == boundaryField_[patchI] + fieldAverage;
///             //         }
///             //     }
///         }
///         else
///         {
///             break;
///         }
///     }
///
/// #endif
///     return ok;
/// }


// processor0/0.004/uniform/lagrangian/kinematicCloud/cloudProperties <stream>   - same for all processors
// processor0/0.004/uniform/lagrangian/kinematicCloud/kinematicCloudOutputProperties <stream>    - same for all processors

// processor0/0.004/lagrangian/kinematicCloud/
// processor0/0.004/lagrangian/kinematicCloud/{U, age, active, angularMom} ...

/// bool Foam::adiosWrite::readScalarFields(adiosReader::helper& helper, label regionID)
/// {
///     bool ok = true;
///
///     const regionInfo& rInfo = regions_[regionID];
///
///     const fieldGroup<scalar>& fields(rInfo.scalarFields_);
///     const fvMesh& mesh = time_.lookupObject<fvMesh>(rInfo.name_);
///
///     forAll(fields, fieldI)
///     {
///         // Lookup field
///         volScalarField& field =
///             const_cast<volScalarField&>
///             (
///                 mesh.lookupObject<volScalarField>(fields[fieldI])
///             );
///
///         Info<< "    readScalarField: " << field.name() << endl;
///
///         fileName datasetName
///         (
///             "region" + Foam::name(regionID)/
///             "fields"/
///             fields[fieldI]
///         );
///
///         {
///             // Read into a plain continuous array for the data
///             ioScalar data[field.size()];
///
///             // Read data from file
///             ok = helper.getDataSet(datasetName, &data);
///             if (ok)
///             {
///                 field.internalField() = UList<ioScalar>(&data[0], field.size());
///             }
///             else
///             {
///                 break;
///             }
///         }
///
///         forAll(field.boundaryField(), patchI)
///         {
///             const fvPatchScalarField& psf = field.boundaryField()[patchI];
///
///             {
///                 Info<< "      patchfield " << patchI
///                     << ": name=" << psf.patch().name()
///                     << ": type=" << psf.type()
///                     << " empty=" << psf.empty()
///                     << " size=" << psf.size()
///                     << endl;
///             }
///
///             if (!isA<emptyFvPatchField<scalar> >(psf))
///             {
///                 continue;
///             }
///
///             ioScalar data[psf.size()];
///
///             ok = helper.getDataSet
///             (
///                 // FIXME: what's the name of psf in the output?
///                 datasetName/"patch" + Foam::name(patchI),
///                 &data
///             );
///
///             // Take reference to scalarField so that we can use = operator
///             // with UList, over-writing fixedValue if present
///             // Note: == operator not available for assignment from UList
///             scalarField& sf = field.boundaryField()[patchI];
///             if (ok)
///             {
///                 sf = UList<ioScalar>(&data[0], sf.size());
///             }
///             else
///             {
///                 break;
///             }
///         }
///
///         if (!ok) break;
///     }
///
///     return ok;
/// }


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

        fileName varPath
        (
            "region" + Foam::name(rInfo.index_)
          / "clouds" / rInfo.cloudNames_[cloudI]
        );
        fileName datasetName(varPath/"nParticlesPerProc");

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

        // Read into a plain continuous array for the data
        // Allocate memory for 1-comp. dataset of type 'integer' for adios_integer reads
        int labelData[nparts];
        // Allocate memory for 1-comp. dataset of type 'float/double' for reads
        //ioScalar scalarData[nparts];

        // Read original processor ID
        if (findStrings(rInfo.cloudAttribs_, "origProc"))
        {
            Info<< "      dataset origProc " << endl;
            ok = helper.getDataSet(varPath/"origProc", labelData);
            label i = 0;
            forAllIter(basicKinematicCloud, *q, pIter)
            {
                pIter().origProc() = labelData[i++];
            }
        }

        if (!ok) break;
    }

    return ok;
}


// ************************************************************************* //
