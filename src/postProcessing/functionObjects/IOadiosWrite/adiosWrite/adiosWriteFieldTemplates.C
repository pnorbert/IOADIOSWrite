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
#include "IOstream.H"
#include "Ostream.H"
#include "OStringStream.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FieldType>
size_t Foam::adiosWrite::fieldDefine
(
    const fvMesh& mesh,
    const regionInfo& rInfo,
    const fieldGroup<typename FieldType::value_type>& fields
)
{
    typedef typename FieldType::value_type pType;

    OStringStream outbuf(IOstream::BINARY);
    size_t maxLen = 0;

    forAll(fields, fieldI)
    {
        // Lookup field
        const FieldType& field = mesh.lookupObject<FieldType>(fields[fieldI]);

        fileName varPath = rInfo.fieldVarPath(fields[fieldI]);

#ifdef FOAM_ADIOS_PATCH_WRITE
        fileName datasetName
        (
            "region" + Foam::name(rInfo.index_)
          / "fields" / fields[fieldI]
        );
#endif
        outbuf.rewind();
        outbuf << field;

        size_t bufLen = outbuf.str().size();
        maxLen = Foam::max(maxLen, bufLen);

        Pout<< "    fieldDefine: " << varPath
            << "  (size " << field.size() << ")"
            << "  stream-size " << bufLen << endl
            << "  stream-content " << outbuf.str() << endl
            << "  ----" << endl;

        adios_define_var
        (
            groupID_,
            varPath.c_str(),                    // name
            "",
            adios_unsigned_byte,                // data-type
            Foam::name(bufLen).c_str(),         // dimensions
            Foam::name(bufLen).c_str(),         // global dimensions
            "0"                                 // local offsets
        );
        outputSize_ += bufLen;

        adios_define_attribute
        (
            groupID_,
            "class",
            varPath.c_str(),
            adios_string,
            field.type().c_str(),  // volScalarField etc.
            NULL
        );

        // could also write field.dimensions() as an attribute
        // if needed to save parsing

#ifdef FOAM_ADIOS_PATCH_WRITE

        // Define a 1D array with field.size as global size, and local (this
        // process') size and with offset 0
        // Type is float or double depending on OpenFoam precision
        adios_define_var
        (
            groupID_,
            datasetName.c_str(),
            "",
            ADIOS_SCALAR,
            Foam::name(field.size()).c_str(),
            Foam::name(field.size()).c_str(),
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
            fileName patchName(datasetName/"patch" + Foam::name(patchI));
            /*
            sprintf(dimstr, "%d", pf1.size()); // global and local dimension of the 1D array
            adios_define_var(groupID_, patchName, "", ADIOS_SCALAR, dimstr, dimstr, "0");
            */
            //adios_define_var(groupID_, patchName, "", adios_string, NULL, NULL, NULL);

            OStringStream ostr;
            ostr<< pf1;
            //Info<< "patch " << patchName << " content: " << ostr.str() << endl;
            adios_define_var
            (
                groupID_,
                patchName.c_str(),
                "",
                adios_byte,
                Foam::name(ostr.str().length()).c_str(),
                Foam::name(ostr.str().length()).c_str(),
                "0"
            );

            // Define attributes to describe this patch
            //char pathstr[128];
            //sprintf(pathstr, "fields/%s/patch%d", fields[fieldI].c_str(), patchI);
            adios_define_attribute
            (
                groupID_,
                "name",
                patchName.c_str(),
                adios_string,
                pf1.patch().name().c_str(),
                NULL
            );
            adios_define_attribute
            (
                groupID_,
                "type",
                patchName.c_str(),
                adios_string,
                pf1.patch().type().c_str(),
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

#endif /* FOAM_ADIOS_PATCH_WRITE */
    }

    Pout<< "max stream-size: " << maxLen << endl;

    return maxLen;
}


template<class FieldType>
void Foam::adiosWrite::fieldWrite
(
    const fvMesh& mesh,
    const regionInfo& rInfo,
    const fieldGroup<typename FieldType::value_type>& fields
)
{
    // typedef typename FieldType::value_type pType;

    OStringStream outbuf(IOstream::BINARY);

    forAll(fields, fieldI)
    {
        // Lookup field
        const FieldType& field = mesh.lookupObject<FieldType>(fields[fieldI]);

        fileName varPath = rInfo.fieldVarPath(fields[fieldI]);

#ifdef FOAM_ADIOS_PATCH_WRITE
        // Dataset for this process
        fileName datasetName
        (
            "region" + Foam::name(rInfo.index_)
          / "fields" / fields[fieldI]
        );
#endif
        outbuf.rewind();
        outbuf << field;

        Info<< "    fieldWrite: " << varPath << endl;

        // Do the actual write
        adios_write
        (
            fileID_,
            varPath.c_str(),
            outbuf.str().c_str()
        );

#ifdef FOAM_ADIOS_PATCH_WRITE

        // ADIOS requires non-const access to the field for writing (?)
        FieldType& fld = const_cast<FieldType&>(field);

        // Do the actual write
        adios_write(fileID_, datasetName.c_str(), fld.data());
        //adios_write (fileID_, datasetName, reinterpret_cast<void *>(fld.internalField().data());

        typename FieldType::GeometricBoundaryField& bf = fld.boundaryField();

        // Do the same for all patches
        forAll(bf, patchI)
        {
            typename FieldType::PatchFieldType& pf1 = bf[patchI];
            Info<< "      patchfield " << patchI << ":" << endl;
#if 0
            // FIXME: what's the name of pf1 in the output?
            sprintf(patchName, "%s/patch%d", datasetName, patchI);
            adios_write(fileID_, patchName, pf1.data());
#else
            OStringStream ostr;
            ostr<< pf1;
            //pf1.write(ostr);
            // Need to copy constant string to buffer because adios_write() wants non-const
            char *buf = strdup(ostr.str().c_str());
            fileName patchName(datasetName/"patch" + Foam::name(patchI));
            Info<< "patch " << patchName << " content: " << ostr.str() << endl;
            adios_write(fileID_, patchName.c_str(), buf);
            free(buf);
#endif
        }

#endif /* FOAM_ADIOS_PATCH_WRITE */

    }
}


// ************************************************************************* //
