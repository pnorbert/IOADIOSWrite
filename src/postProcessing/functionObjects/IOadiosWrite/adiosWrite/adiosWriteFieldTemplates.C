/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |               2015 Norbert Podhorszki
                            |               2016 OpenCFD Ltd.
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
#include "OStringStream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FieldType>
void Foam::adiosWrite::fieldDefine
(
    const fvMesh& mesh,
    const fieldGroup<typename FieldType::value_type>& fields,
    const label regionID
)
{
    typedef typename FieldType::value_type pType;

    char datasetName[256];
    char patchName[256];
    char dimstr[16];
    forAll(fields, fieldI)
    {
        Info<< "    fieldDefine: " << fields[fieldI] << endl;

        // Lookup field
        const FieldType& field = mesh.lookupObject<FieldType>(fields[fieldI]);

        /*sprintf
            (
                datasetName,
                "FIELDS/%s/processor%i/%s",
                mesh_.time().timeName().c_str(),
                Pstream::myProcNo(),
                fields[fieldI].c_str()
            );*/
        sprintf
        (
            datasetName,
            "region%d/fields/%s",
            regionID,
            fields[fieldI].c_str()
        );

        // Define a 1D array with field.size as global size, and local (this
        // process') size and with offset 0
        // Type is float or double depending on OpenFoam precision
        sprintf(dimstr, "%d", field.size()); // global and local dimension of the 1D array
        adios_define_var
        (
            groupID_,
            datasetName,
            "",
            ADIOS_SCALAR,
            dimstr,
            dimstr,
            "0"
        );

        // Count the total size we are going to write from this process
        outputSize_ += field.size()*pTraits<pType>::nComponents*sizeof(ioScalar);

        const typename FieldType::GeometricBoundaryField& bf =
            field.boundaryField();

        forAll(bf, patchI)
        {
            const typename FieldType::PatchFieldType& pf1 = bf[patchI];
            Info<< "      patchfield " << patchI
                << ": name=" << pf1.patch().name()
                << ": type=" << pf1.type()
                << " empty=" << pf1.empty()
                << " size=" << pf1.size()
                << endl;
            // FIXME: what's the name of pf1 in the output?
            sprintf(patchName, "%s/patch%d", datasetName, patchI);
            /*
            sprintf(dimstr, "%d", pf1.size()); // global and local dimension of the 1D array
            adios_define_var(groupID_, patchName, "", ADIOS_SCALAR, dimstr, dimstr, "0");
            */
            //adios_define_var(groupID_, patchName, "", adios_string, NULL, NULL, NULL);

            OStringStream ostr;
            ostr<< pf1;
            //Info<< "patch " << patchName << " content: " << ostr.str() << endl;
            sprintf(dimstr, "%d", ostr.str().length()); // global and local dimension of the 1D array
            adios_define_var
            (
                groupID_,
                patchName,
                "",
                adios_byte,
                dimstr,
                dimstr,
                "0"
            );

            // Define attributes to describe this patch
            //char pathstr[128];
            char tmpstr[128];
            //sprintf(pathstr, "fields/%s/patch%d", fields[fieldI].c_str(), patchI);
            sprintf(tmpstr, "%s", pf1.patch().name().c_str());
            adios_define_attribute
            (
                groupID_,
                "name",
                patchName,
                adios_string,
                tmpstr,
                NULL
            );
            sprintf(tmpstr, "%s", pf1.type().c_str());
            adios_define_attribute
            (
                groupID_,
                "type",
                patchName,
                adios_string,
                tmpstr,
                NULL
            );

            // Count the total size we are going to write from this process
            //outputSize_ += pf1.size()*pType::nComponents*sizeof(ioScalar); // size of data array
            outputSize_ += ostr.str().length(); // size of string storage of patch, unknown

            //forAll(pf1, faceI)
            // {
            //     pf1[faceI] = ...
            // }
        }
    }
}


template<class FieldType>
void Foam::adiosWrite::fieldWrite
(
    const fvMesh& mesh,
    const fieldGroup<typename FieldType::value_type>& fields,
    const label regionID
)
{
    typedef typename FieldType::value_type pType;

    forAll(fields, fieldI)
    {
        Info<< "    fieldWriteScalar: " << fields[fieldI] << endl;

        // Lookup field
        const FieldType& field = mesh.lookupObject<FieldType>(fields[fieldI]);

        // Initialize a plain continuous array for the data
        ioScalar *scalarData;
        // This does not work since adios_write() cannot accept const void *
        //scalarData = reinterpret_cast<ioScalar*>(field.internalField().cdata());

        /* FIXME: if we can get the field.xxxx as the double/float array
         * then there is no need for this element-by-element copy
         */
        scalarData = new ioScalar[field.size()*pTraits<pType>::nComponents];
        //if (contiguous("Type of field which is a templated class' object"))  //
        //{
            // adios_write() does not accept const void * arrays, so we need to copy the cdata here
            // This only works as long as the field's type is the same as ioScalar
            // i.e. float or double depending on the WM_DP compile flag
            memcpy
            (
                scalarData,
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

        // Dataset for this process
        sprintf
        (
            datasetName,
            "region%d/fields/%s",
            regionID,
            fields[fieldI].c_str()
        );

        // Do the actual write
        adios_write(fileID_, datasetName, scalarData);
        //adios_write (fileID_, datasetName, reinterpret_cast<void *>(field.internalField().data());

        // Release memory
        delete [] scalarData;

        const typename FieldType::GeometricBoundaryField& bf =
            field.boundaryField();

        // Do the same for all patches
        forAll(bf, patchI)
        {
            const typename FieldType::PatchFieldType& pf1 = bf[patchI];
            Info<< "      patchfield " << patchI << ":" << endl;
#if 0
            scalarData = new ioScalar[field.size()];
            /*memcpy (scalarData,
                    reinterpret_cast<const char *>(pf1.internalField().cdata()),
                    pf1.byteSize()
                   );*/
            forAll(pf1, iter)
            {
                scalarData[iter] = pf1[iter];
            }

            // FIXME: what's the name of pf1 in the output?
            sprintf (patchName, "%s/patch%d", datasetName, patchI);
            adios_write (fileID_, patchName, scalarData);
            delete [] scalarData;
#else
            OStringStream ostr;
            ostr<< pf1;
            //pf1.write(ostr);
            // Need to copy constant string to buffer because adios_write() wants non-const
            char *buf = strdup(ostr.str().c_str());
            sprintf(patchName, "%s/patch%d", datasetName, patchI);
            Info<< "patch " << patchName << " content: " << ostr.str() << endl;
            adios_write(fileID_, patchName, buf);
            free(buf);
#endif
        }
    }
}


// ************************************************************************* //
