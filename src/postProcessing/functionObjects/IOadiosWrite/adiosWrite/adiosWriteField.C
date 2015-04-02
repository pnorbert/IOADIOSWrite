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

void Foam::adiosWrite::fieldDefine()
{
    Info<< "  adiosWrite::fieldDefine:"  << endl;
    fieldDefineScalar();
    fieldDefineVector();
}

void Foam::adiosWrite::fieldWrite()
{
    Info<< "  adiosWrite::fieldWrite:"  << endl;
    fieldWriteScalar();
    fieldWriteVector();
}

void Foam::adiosWrite::fieldDefineScalar()
{
    char datasetName[128];
    char dimstr[16];
    forAll(scalarFields_, fieldI)
    {
        Info<< "    fieldDefineScalar: " << scalarFields_[fieldI] << endl;

        // Lookup field
        const volScalarField& field = obr_.lookupObject<volScalarField>
            (
                scalarFields_[fieldI]
            );

        /*sprintf
            (
                datasetName,
                "FIELDS/%s/processor%i/%s",
                mesh_.time().timeName().c_str(),
                Pstream::myProcNo(),
                scalarFields_[fieldI].c_str()
            );*/
        sprintf (datasetName, "fields/%s", scalarFields_[fieldI].c_str());

        // Define a 1D array with field.size as global size, and 
        //   local (this process') size and with offset 0, 
        // Type is float or double depending on OpenFoam precision
        sprintf (dimstr, "%d", field.size()); // global and local dimension of the 1D array
        adios_define_var (groupID_, datasetName, "", ADIOS_SCALAR, dimstr, dimstr, "0");

        // count the total size we are going to write from this process
        outputSize_ += field.size() * sizeof(ioScalar);
    }
}

        
void Foam::adiosWrite::fieldDefineVector()
{
    char datasetName[128];
    char dimstr[16];
    forAll(vectorFields_, fieldI)
    {
        Info<< "    fieldDefineVector: " << vectorFields_[fieldI] << endl;

        // Lookup field
        const volVectorField& field = obr_.lookupObject<volVectorField>
            (
                vectorFields_[fieldI]
            );

        /*sprintf
            (
                datasetName,
                "FIELDS/%s/processor%i/%s",
                mesh_.time().timeName().c_str(),
                Pstream::myProcNo(),
                vectorFields_[fieldI].c_str()
            );*/
        sprintf (datasetName, "fields/%s", vectorFields_[fieldI].c_str());

        // Define a 2D array with "field.size x 3"  as global size, and 
        //   local (this process') size and with offset 0,0 
        // Type is float or double depending on OpenFoam precision
        sprintf (dimstr, "%d,3", field.size()); // global and local dimension of the 1D array
        adios_define_var (groupID_, datasetName, "", ADIOS_SCALAR, dimstr, dimstr, "0,0");

        // count the total size we are going to write from this process
        outputSize_ += field.size() * 3 * sizeof(ioScalar);
    }
}

void Foam::adiosWrite::fieldWriteScalar()
{
    forAll(scalarFields_, fieldI)
    {
        Info<< "    fieldWriteScalar: " << scalarFields_[fieldI] << endl;
        
        // Lookup field
        const volScalarField& field = obr_.lookupObject<volScalarField>
            (
                scalarFields_[fieldI]
            );
        
        
        // Initialize a plain continuous array for the data
        ioScalar* scalarData;
        scalarData = new ioScalar[field.size()];
        /* FIXME: if we can get the field.xxxx as the double/float array
         * then there is no need for this element-by-element copy
         */
        //if (field.contiguous()) 
        //{
           memcpy (scalarData,  
                   reinterpret_cast<const char *>(field.cdata()),
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

        char datasetName[80];
        // dataset for this process
        /*sprintf
            (
                datasetName,
                "FIELDS/%s/processor%i/%s",
                mesh_.time().timeName().c_str(),
                Pstream::myProcNo(),
                scalarFields_[fieldI].c_str()
            );*/
        sprintf (datasetName, "fields/%s", scalarFields_[fieldI].c_str());

        // Do the actual write
        adios_write (fileID_, datasetName, scalarData);

        // Release memory
        delete [] scalarData;
    }
}


void Foam::adiosWrite::fieldWriteVector()
{
    forAll(vectorFields_, fieldI)
    {
        Info<< "    fieldWriteVector: " << vectorFields_[fieldI] << endl;
        
        const volVectorField& field = obr_.lookupObject<volVectorField>
            (
                vectorFields_[fieldI]
            );
        
        //Initialize a plain continous array for the data
        ioScalar* vectorData;
        vectorData = new ioScalar[field.size()*3];
        
        
        // Loop through the field and construct the array
        forAll(field, iter)
        {
            vectorData[3*iter+0] = field[iter].x();
            vectorData[3*iter+1] = field[iter].y();
            vectorData[3*iter+2] = field[iter].z();
        }
        
        
        // Create the different datasets (needs to be done collectively)
        char datasetName[80];

        // Dataset for this process
        /*sprintf
            (
                datasetName,
                "FIELDS/%s/processor%i/%s",
                mesh_.time().timeName().c_str(),
                Pstream::myProcNo(),
                vectorFields_[fieldI].c_str()
            );*/
        sprintf (datasetName, "fields/%s", vectorFields_[fieldI].c_str());
        
        // Do the actual write
        adios_write (fileID_, datasetName, vectorData);
        
        
        // Release memory
        delete [] vectorData;
    }
}


// ************************************************************************* //
