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
#include "dictionary.H"
#include "scalar.H"
#include "basicKinematicCloud.H"
#include "emptyFvPatchField.H"
#include "IStringStream.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

static ADIOS_FILE *f;

bool Foam::adiosWrite::readData(const fileName& bpFile) const
{
    Info<< " Read data of step " << bpFile << endl;
    f = adios_read_open_file
    (
        bpFile.c_str(),
        ADIOS_READ_METHOD_BP,
        comm_
    );

    if (!f)
    {
        return false;
    }

    Info<< "found num-vars: " << f->nvars << endl;

    size_t maxLen = 0;

    // TODO? restrict sizing to current processor!
    for (int varI=0; varI < f->nvars; ++varI)
    {
        const char * varName = f->var_namelist[varI];
        ADIOS_VARINFO *varInfo = adios_inq_var(f, varName);

        if (!varInfo)
        {
            WarningInFunction
                << "Error reading variable information " << varName
                << " from adios file: "
                << adios_errmsg() << endl;
            continue;
        }

        size_t bytes = 1; // fallback to scalar
        if (varInfo->ndim > 0)
        {
            // only consider 1-D storage:
            bytes = varInfo->dims[0];
            // for (int dimI=1; dimI < varInfo->ndim; ++dimI)
            // {
            //     ... varInfo->dims[dimI];
            // }
        }

        bytes *= adios_type_size(varInfo->type, const_cast<char *>(""));

        Info<< "   >" << varName << endl;
        Info<< "   >>" << adios_type_to_string(varInfo->type) << endl;
        Info<< "   bytes>" << bytes << endl;

        maxLen = Foam::max(maxLen, bytes);

        // free ADIOS_VARINFO
        adios_free_varinfo(varInfo);
    }


    Info<< "   max-length: " << maxLen << endl;

    IStringStream is("", IOstream::BINARY);    // empty buffer to avoid copy-in behaviour
    char chardata[maxLen];                     // this is our real data buffer

    // attach our buffer
    is.stdStream().rdbuf()->pubsetbuf(chardata, maxLen);

    // Process N will read writeblock N
    ADIOS_SELECTION *sel = adios_selection_writeblock(Pstream::myProcNo());
    {

        fileName varName
        (
            "region" + Foam::name(0)
          / "field" / "pMean" / "stream"
        );

        // Read into a plain continuous array for the data
        // char data[maxLen];

        Info<<"read back in " << varName << endl;

        is.rewind();
        int err = adios_schedule_read(f, sel, varName.c_str(), 0, 1, chardata);

        if (err)
        {
            WarningInFunction
                << "Error reading field " << varName
                << " from adios file: "
                << adios_errmsg() << endl;
        }
        else
        {
            err = adios_perform_reads(f, 1); // blocking
        }

        if (!err)
        {
            Info<<"OK" << endl;

            char c;
            for (int i=0; i < 100; ++i)
            {
                is.get(c);
                Info<< "char: " << c << endl;

                chardata[i] = c ^ '#';
            }

            is.rewind();
            for (int i=0; i < 100; ++i)
            {
                is.get(c);
                Info<< "updated: " << c << endl;

                chardata[i] = c ^ '#';
            }

            is.rewind();
            for (int i=0; i < 100; ++i)
            {
                is.get(c);
                Info<< "again: " << c << endl;
            }
        }
    }

    adios_selection_delete(sel);

    // fileNameList varList(f->nvars);

    bool ok = true;
    forAll(regions_, regionID)
    {
        ok =
        (
            ok
         && readScalarFields(regionID)
         // && readSurfaceScalarFields(regionID)
         // && readVectorFields(regionID);
         && readClouds(regionID)
        );
    }

    adios_read_close(f);
    return ok;
}


bool Foam::adiosWrite::readData(const instant& when) const
{
    Info<< " Read data of step " << when.name() << endl;
    return readData(when.name());
}


bool Foam::adiosWrite::readData() const
{
    Info<< " Read data of step " << obr_.time().timeName() << endl;
    return readData(dataDirectory/obr_.time().timeName() + ".bp");
}


bool Foam::adiosWrite::readADIOSVar
(
    ADIOS_FILE *f,
    ADIOS_SELECTION *sel,
    const char *path1,
    const char *path2,
    void *data
) const
{
    fileName datasetName(path1);

    if (path2)
    {
        datasetName = datasetName/path2;
    }

    // IstreamAdios(f, sel, name);
    bool ok = true;
    int err = adios_schedule_read(f, sel, datasetName.c_str(), 0, 1, data);
    if (err)
    {
        WarningInFunction
            << "Error reading field " << datasetName
            << " from adios checkpoint file: "
            << adios_errmsg()
            << endl;
        ok = false;
    }
    else
    {
        err = adios_perform_reads(f, 1); // blocking
    }
    return ok;
}


bool Foam::adiosWrite::readScalarFields(label regionID) const
{
    bool succ = true;

    // Process N will read writeblock N
    ADIOS_SELECTION *sel = adios_selection_writeblock(Pstream::myProcNo());

    const fieldGroup<scalar>& fields(regions_[regionID].scalarFields_);
    const fvMesh& mesh = time_.lookupObject<fvMesh>(regions_[regionID].name_);
    forAll(fields, fieldI)
    {
        // Lookup field
        volScalarField& field =
            const_cast<volScalarField&>
            (
                mesh.lookupObject<volScalarField>(fields[fieldI])
            );

        Info<< "    readScalarField: " << field.name() << endl;


        fileName varName
        (
            "region" + Foam::name(regionID)
          / "field" / fields[fieldI] / "stream"
        );


        fileName datasetName
        (
            "region" + Foam::name(regionID)/
            "fields"/
            fields[fieldI]
        );

        {
            // Read into a plain continuous array for the data
            ioScalar data[field.size()];

            // Read data from file
            succ = readADIOSVar(f, sel, datasetName.c_str(), NULL, &data);
            if (succ)
            {
                field.internalField() = UList<ioScalar>(&data[0], field.size());
            }
            else
            {
                break;
            }
        }

        forAll(field.boundaryField(), patchI)
        {
            const fvPatchScalarField& psf = field.boundaryField()[patchI];

            {
                Info<< "      patchfield " << patchI
                    << ": name=" << psf.patch().name()
                    << ": type=" << psf.type()
                    << " empty=" << psf.empty()
                    << " size=" << psf.size()
                    << endl;
            }

            if (!isA<emptyFvPatchField<scalar> >(psf))
            {
                continue;
            }

            // FIXME: what's the name of psf in the output?
            fileName patchName(datasetName/"patch" + Foam::name(patchI));

            ioScalar data[psf.size()];

            succ = readADIOSVar(f, sel, patchName.c_str(), NULL, &data);
            // Take reference to scalarField so that we can use = operator
            // with UList, over-writing fixedValue if present
            // Note: == operator not available for assignment from UList
            scalarField& sf = field.boundaryField()[patchI];
            if (succ)
            {
                sf = UList<ioScalar>(&data[0], sf.size());
            }
            else
            {
                break;
            }
        }

        if (!succ) break;
    }

    adios_selection_delete(sel);
    return succ;
}



bool Foam::adiosWrite::readClouds(label regionID) const
{
    bool succ = true;

    // Process N will read writeblock N
    ADIOS_SELECTION *sel = adios_selection_writeblock(Pstream::myProcNo());

    const regionInfo& rInfo = regions_[regionID];
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
            "region" + Foam::name(regionID)/
            "clouds"/
            rInfo.cloudNames_[cloudI]
        );
        fileName datasetName(varPath/"nParticlesPerProc");

        // Get the number of particles saved by this rank
        int nparts = 0;
        succ = readADIOSVar(f, sel, datasetName.c_str(), NULL, &nparts);

        /*
        // Get the number of particles saved by this rank
        ADIOS_VARINFO * vi = adios_inq_var(f, datasetName);
        if (vi != NULL) {
            adios_inq_var_blockinfo(f, vi);
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
            succ = false;
        }
        */

        if (!succ) break;

        // Read into a plain continuous array for the data
        // Allocate memory for 1-comp. dataset of type 'integer' for adios_integer reads
        int labelData[nparts];
        // Allocate memory for 1-comp. dataset of type 'float/double' for reads
        //ioScalar scalarData[nparts];

        // Read original processor ID
        if (findStrings(rInfo.cloudAttribs_, "origProc"))
        {
            Info<< "      dataset origProc " << endl;
            //sprintf(datasetName, "%s/origProc", varPath);
            succ = readADIOSVar(f, sel, varPath.c_str(), "origProc", labelData);
            label i = 0;
            forAllIter(basicKinematicCloud, *q, pIter)
            {
                pIter().origProc() = labelData[i++];
            }
        }

        if (!succ) break;
    }

    adios_selection_delete(sel);
    return succ;
}


// ************************************************************************* //
