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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::adiosWrite::fieldDefine(label regionID)
{
    char tmpstr[128];
    Info<< "  adiosWrite::fieldDefine: region " << regionID << " " << regions_[regionID].name_ << ": " << endl;
    sprintf (tmpstr, "region%d", regionID);
    adios_define_attribute (groupID_, "name", tmpstr, adios_string, regions_[regionID].name_.c_str(), NULL);

    fieldDefineScalar(regionID);
    fieldDefineVector(regionID);
}

void Foam::adiosWrite::fieldWrite(label regionID)
{
    Info<< "  adiosWrite::fieldWrite: region " << regionID << " " << regions_[regionID].name_ << ": " << endl;
    fieldWriteScalar(regionID);
    fieldWriteVector(regionID);
}

void Foam::adiosWrite::fieldDefineScalar(label regionID)
{
    char datasetName[256];
    char patchName[256];
    char dimstr[16];
    const fieldGroup<scalar>& fields (regions_[regionID].scalarFields_);
    const fvMesh& m = time_.lookupObject<fvMesh>(regions_[regionID].name_);
    forAll(fields, fieldI)
    {
        Info<< "    fieldDefineScalar: " << fields[fieldI] << endl;

        // Lookup field
        const volScalarField& field = m.lookupObject<volScalarField>
            (
                fields[fieldI]
            );

        /*sprintf
            (
                datasetName,
                "FIELDS/%s/processor%i/%s",
                mesh_.time().timeName().c_str(),
                Pstream::myProcNo(),
                fields[fieldI].c_str()
            );*/
        sprintf (datasetName, "region%d/fields/%s", regionID, fields[fieldI].c_str());

        // Define a 1D array with field.size as global size, and 
        //   local (this process') size and with offset 0, 
        // Type is float or double depending on OpenFoam precision
        sprintf (dimstr, "%d", field.size()); // global and local dimension of the 1D array
        adios_define_var (groupID_, datasetName, "", ADIOS_SCALAR, dimstr, dimstr, "0");

        // count the total size we are going to write from this process
        outputSize_ += field.size() * sizeof(ioScalar);

        forAll(field.boundaryField(), patchI)
        {
            const fvPatchScalarField& pf1 = field.boundaryField()[patchI];
            Info<< "      patchfield " << patchI 
                << ": name=" << pf1.patch().name() 
                << ": type=" << pf1.type() 
                << " empty=" << pf1.empty() 
                << " size=" << pf1.size() 
                << endl;
            // FIXME: what's the name of pf1 in the output?
            sprintf (patchName, "%s/patch%d", datasetName, patchI);
            sprintf (dimstr, "%d", pf1.size()); // global and local dimension of the 1D array
            adios_define_var (groupID_, patchName, "", ADIOS_SCALAR, dimstr, dimstr, "0");

            // define attributes to describe this patch
            //char pathstr[128];
            char tmpstr[128];
            //sprintf (pathstr, "fields/%s/patch%d", fields[fieldI].c_str(), patchI);
            sprintf (tmpstr, "%s", pf1.patch().name().c_str());
            adios_define_attribute (groupID_, "name", patchName, adios_string, tmpstr, NULL);
            sprintf (tmpstr, "%s", pf1.type().c_str());
            adios_define_attribute (groupID_, "type", patchName, adios_string, tmpstr, NULL);

            // count the total size we are going to write from this process
            outputSize_ += pf1.size() * sizeof(ioScalar);

            //forAll(pf1, faceI)
            // {
            //     pf1[faceI] = ...
            // }
        }

    }
}

        
void Foam::adiosWrite::fieldDefineVector(label regionID)
{
    char datasetName[256];
    char dimstr[16];
    const fieldGroup<vector>& fields (regions_[regionID].vectorFields_);
    const fvMesh& m = time_.lookupObject<fvMesh>(regions_[regionID].name_);
    forAll(fields, fieldI)
    {
        Info<< "    fieldDefineVector: " << fields[fieldI] << endl;

        // Lookup field
        const volVectorField& field = m.lookupObject<volVectorField>
            (
                fields[fieldI]
            );

        /*sprintf
            (
                datasetName,
                "FIELDS/%s/processor%i/%s",
                mesh_.time().timeName().c_str(),
                Pstream::myProcNo(),
                fields[fieldI].c_str()
            );*/
        sprintf (datasetName, "region%d/fields/%s", regionID, fields[fieldI].c_str());

        // Define a 2D array with "field.size x 3"  as global size, and 
        //   local (this process') size and with offset 0,0 
        // Type is float or double depending on OpenFoam precision
        sprintf (dimstr, "%d,3", field.size()); // global and local dimension of the 1D array
        adios_define_var (groupID_, datasetName, "", ADIOS_SCALAR, dimstr, dimstr, "0,0");

        // count the total size we are going to write from this process
        outputSize_ += field.size() * 3 * sizeof(ioScalar);
    }
}

void Foam::adiosWrite::fieldWriteScalar(label regionID)
{
    const fieldGroup<scalar>& fields (regions_[regionID].scalarFields_);
    const fvMesh& m = time_.lookupObject<fvMesh>(regions_[regionID].name_);
    forAll(fields, fieldI)
    {
        Info<< "    fieldWriteScalar: " << fields[fieldI] << endl;
        
        // Lookup field
        const volScalarField& field = m.lookupObject<volScalarField>
            (
                fields[fieldI]
            );
        
        
        // Initialize a plain continuous array for the data
        ioScalar* scalarData;
        // This does not work since adios_write() cannot accept const void *
        //scalarData = reinterpret_cast<ioScalar*>(field.internalField().cdata());

        /* FIXME: if we can get the field.xxxx as the double/float array
         * then there is no need for this element-by-element copy
         */
        scalarData = new ioScalar[field.size()];
        //if (contiguous("Type of field which is a templated class' object"))  //
        //{
            // adios_write() does not accept const void * arrays, so we need to copy the cdata here 
            // This only works as long as the field's type is the same as ioScalar
            // i.e. float or double depending on the WM_DP compile flag
           memcpy (scalarData,  
                   //reinterpret_cast<const char *>(field.cdata()),
                   reinterpret_cast<const char *>(field.internalField().cdata()),
                   field.byteSize()
                  );
        //}
        //else
        //{
        //    // Loop through the field and construct the array
        //    forAll(field, iter)
        //    {
        //        scalarData[iter] = field[iter];
        //    }
        //}

        char datasetName[256];
        char patchName[256];
        // dataset for this process
        sprintf (datasetName, "region%d/fields/%s", regionID, fields[fieldI].c_str());

        // Do the actual write
        adios_write (fileID_, datasetName, scalarData);
        //adios_write (fileID_, datasetName, reinterpret_cast<void *>(field.internalField().data());

        // Release memory
        delete [] scalarData;


        // Do the same for all patches
        forAll(field.boundaryField(), patchI)
        {
            const fvPatchScalarField& pf1 = field.boundaryField()[patchI];
            Info<< "      patchfield " << patchI << ":" << endl;
            scalarData = new ioScalar[field.size()];
            memcpy (scalarData,  
                    reinterpret_cast<const char *>(pf1.internalField().cdata()),
                    pf1.byteSize()
                   );
            // FIXME: what's the name of pf1 in the output?
            sprintf (patchName, "%s/patch%d", datasetName, patchI);
            adios_write (fileID_, patchName, scalarData);
            delete [] scalarData;
        }
    }
}


void Foam::adiosWrite::fieldWriteVector(label regionID)
{
    const fieldGroup<vector>& fields (regions_[regionID].vectorFields_);
    const fvMesh& m = time_.lookupObject<fvMesh>(regions_[regionID].name_);
    forAll(fields, fieldI)
    {
        Info<< "    fieldWriteVector: " << fields[fieldI] << endl;
        
        const volVectorField& field = m.lookupObject<volVectorField>
            (
                fields[fieldI]
            );
        
        //Initialize a plain continous array for the data
        ioScalar* vectorData;
        vectorData = new ioScalar[field.size()*3];
        
        //Info<< "      vector size = " << field.size() 
        //    << "  byte size = " << field.byteSize() << endl;
        
        //if (field.contiguous()) 
        //{
           memcpy (vectorData,  
                   //reinterpret_cast<const char *>(field.cdata()),
                   reinterpret_cast<const char *>(field.internalField().cdata()),
                   field.byteSize()
                  );
        //}
        //else
        //{
        // Loop through the field and construct the array
        //forAll(field, iter)
        //{
        //    vectorData[3*iter+0] = field[iter].x();
        //    vectorData[3*iter+1] = field[iter].y();
        //    vectorData[3*iter+2] = field[iter].z();
        //}
        //}
        
        // Create the different datasets (needs to be done collectively)
        char datasetName[256];
        char patchName[256];

        // Dataset for this process
        sprintf (datasetName, "region%d/fields/%s", regionID, fields[fieldI].c_str());
        
        // Do the actual write
        adios_write (fileID_, datasetName, vectorData);
        
        
        // Release memory
        delete [] vectorData;
    }
}


// ************************************************************************* //
