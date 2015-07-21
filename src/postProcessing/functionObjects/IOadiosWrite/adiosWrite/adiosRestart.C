/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |               2015 Norbert Podhorszki
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



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

static ADIOS_FILE * f;

bool Foam::adiosWrite::readData(scalar time)
{
    Info<< " Read data of step " << time << endl;
    bool succ = true;
    char fname[256];
    snprintf (fname, sizeof(fname), "adiosData/%g.bp", time);
    f = adios_read_open_file (fname, ADIOS_READ_METHOD_BP, comm_);
    if (!f) {
        return false;
    }

    forAll(regions_, regionID)
    {
        if (succ) {
            succ = readScalarFields(regionID);
        }
        if (succ) {
            //succ = readSurfaceScalarFields(regionID);
        }
        if (succ) {
            //succ = readVectorFields(regionID);
        }
        if (succ) {
            succ = readClouds(regionID);
        }
    }

    adios_read_close (f);
    return succ;
}

bool Foam::adiosWrite::readADIOSVar (
        ADIOS_FILE *f, ADIOS_SELECTION *sel, 
        const char * path1, const char * path2, 
        void * data)
{
    char datasetName[256];
    int err;
    bool succ = true;
    if (path2) 
    {
        sprintf (datasetName, "%s/%s", path1, path2);
    } 
    else 
    {
        sprintf (datasetName, "%s", path1);
    }

    err = adios_schedule_read (f, sel, datasetName, 0, 1, data);
    if (!err) {
        err = adios_perform_reads (f, 1);
    }
    if (err) {
        WarningIn ("Foam::adiosWrite::read_from_adios()")  
            << "Error reading field " << datasetName 
            << " from adios checkpoint file: " 
            << adios_errmsg() 
            << endl;
        succ = false;
    }
    return succ;
}


bool Foam::adiosWrite::readScalarFields(label regionID)
{
    bool succ = true;
    char datasetName[256];
    char patchName[256];

    // process N will read writeblock N
    ADIOS_SELECTION * sel = adios_selection_writeblock (Pstream::myProcNo());

    const fieldGroup<scalar>& fields (regions_[regionID].scalarFields_);
    const fvMesh& m = time_.lookupObject<fvMesh>(regions_[regionID].name_);
    forAll(fields, fieldI)
    {
        Info<< "    readScalarField: " << fields[fieldI] << endl;

        // Lookup field
        const volScalarField& field = m.lookupObject<volScalarField>
            (
             fields[fieldI]
            );

        sprintf (datasetName, "region%d/fields/%s", regionID, fields[fieldI].c_str());

        // Read into a plain continuous array for the data
        ioScalar* scalarData;
        scalarData = new ioScalar[field.size()];

        // read data from file
        succ = readADIOSVar (f, sel, datasetName, NULL, scalarData);
        if (succ) 
        {
            /* FIXME: Need to copy the data somehow into
               field.internalField().data()
             */
        }

        // Release temporary memory
        delete [] scalarData;

        if (!succ) break;

        forAll(field.boundaryField(), patchI)
        {
            char buffer[65536];
            const fvPatchScalarField& pf1 = field.boundaryField()[patchI];
            Info<< "      patchfield " << patchI
                << ": name=" << pf1.patch().name()
                << ": type=" << pf1.type()
                << " empty=" << pf1.empty()
                << " size=" << pf1.size()
                << endl;
            // FIXME: what's the name of pf1 in the output?
            sprintf (patchName, "%s/patch%d", datasetName, patchI);

            succ = readADIOSVar (f, sel, patchName, NULL, scalarData);

            if (!succ) break;

            IStringStream istr (buffer);
            /* FIXME: Need to update the patch with this information now */
            /*
            pf1 << istr;
            */

            Info<< "      patchfield " << patchI
                << ": name=" << pf1.patch().name()
                << ": type=" << pf1.type()
                << " empty=" << pf1.empty()
                << " size=" << pf1.size()
                << " after reading " << istr.str().length() << " bytes"
                << endl;
        }

        if (!succ) break;
    }

    adios_selection_delete (sel);
    return succ;
}



bool Foam::adiosWrite::readClouds(label regionID)
{
    bool succ = true;
    char varPath[256];
    char datasetName[256];
    char patchName[256];

    // process N will read writeblock N
    ADIOS_SELECTION * sel = adios_selection_writeblock (Pstream::myProcNo());

    regionInfo& r = regions_[regionID];
    const fvMesh& m = time_.lookupObject<fvMesh>(regions_[regionID].name_);
    forAll(r.cloudNames_, cloudI)
    {
        Info<< "    cloud: " << r.cloudNames_[cloudI] << endl;

        const kinematicCloud& cloud = m.lookupObject<kinematicCloud>
            (
             r.cloudNames_[cloudI]
            );
        //basicKinematicCloud *q = (basicKinematicCloud*) &cloud;
        const basicKinematicCloud *q = reinterpret_cast<const basicKinematicCloud*>(&cloud);


        sprintf (varPath, "region%d/clouds/%s", regionID, r.cloudNames_[cloudI].c_str());

        // get the number of particles saved by this rank
        int  nparts = 0;
        sprintf (datasetName, "%s/nParticlesPerProc", varPath);
        succ = readADIOSVar (f, sel, datasetName, NULL, &nparts);

        /*
        // get the number of particles saved by this rank
        ADIOS_VARINFO * vi = adios_inq_var (f, datasetName);
        if (vi != NULL) {
            adios_inq_var_blockinfo (f, vi);
            if (vi->sum_nblocks > Pstream::myProcNo())
            {
                nParticles = vi->blockinfo[Pstream::myProcNo()].count[0];
            }
            else
            {
                WarningIn ("Foam::adiosWrite::readClouds()")  
                    << "Error reading cloud " << datasetName 
                    << ": Number of available blocks = " << vi->sum_nblocks
                    << " is less than this rank " << Pstream::myProcNo() 
                    << endl;
            }
            adios_free_varinfo (vi);
        } 
        else 
        {
            WarningIn ("Foam::adiosWrite::readClouds()")  
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
        ioScalar* scalarData;
        scalarData = new ioScalar[nparts];

        // Read original processor ID
        if (findStrings(r.cloudAttribs_, "origProc"))
        {
            Info<< "      dataset origProc " << endl;
            //sprintf (datasetName, "%s/origProc", varPath);
            succ = readADIOSVar (f, sel, varPath, "origProc", labelData);
            /* FIXME: Need to copy the data somehow into
               label i = 0;
               forAllIter(basicKinematicCloud, *q, pIter)
               {
                   pIter().origProc() <--- labelData[i];
                   i++;
               }
             */
        }
        if (!succ) break;


        // Release temporary memory
        delete [] scalarData;


        if (!succ) break;
    }

    adios_selection_delete (sel);
    return succ;
}





