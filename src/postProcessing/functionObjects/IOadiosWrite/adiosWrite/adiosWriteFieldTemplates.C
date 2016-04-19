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

#include "OCountStream.H"
#include "OCompactCountStream.H"
#include "OCompactStringStream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FieldType>
size_t Foam::adiosWrite::fieldDefine
(
    const fvMesh& mesh,
    const regionInfo& rInfo,
    const fieldGroup<typename FieldType::value_type>& fields
)
{
#ifdef FOAM_ADIOS_PATCH_WRITE
    typedef typename FieldType::value_type pType;
#endif
    // OCountStream os(adiosCore::strFormat);
    OCompactCountStream os(adiosCore::strFormat);

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
        os.rewind();
        os << field;
        size_t bufLen = os.size();
        maxLen = Foam::max(maxLen, bufLen);

        // OStringStream check(adiosCore::strFormat);
        OCompactStringStream check(adiosCore::strFormat);
        check << field;
        size_t oslen = check.str().size();

        Pout<< "    fieldDefine: " << varPath
            << "  (size " << field.size() << ")"
            << "  stream-size " << oslen << " (counted " << bufLen << ")" << endl
            << "  stream-content " << check.str().c_str() << endl
            << "  ----" << endl;

        defineVariable(varPath, adios_unsigned_byte, bufLen);

        // volScalarField etc.
        defineAttribute("class", varPath, field.type());

        // could also write field.dimensions() as an attribute
        // if needed to save parsing

#ifdef FOAM_ADIOS_PATCH_WRITE

        // Define a 1D array with field.size as global size, and local (this
        // process') size and with offset 0
        // Type is float or double depending on OpenFoam precision
        defineVariable(datasetName, ADIOS_SCALAR, field.size());

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

            OStringStream ostr(adiosCore::strFormat);
            ostr<< pf1;
            //Info<< "patch " << patchName << " content: " << ostr.str() << endl;
            defineVariable(patchName, adios_unsigned_byte, ostr.str().length());

            // Define attributes to describe this patch
            //char pathstr[128];
            //sprintf(pathstr, "fields/%s/patch%d", fields[fieldI].c_str(), patchI);
            defineAttribute("name", patchName, pf1.patch().name());
            defineAttribute("type", patchName, pf1.patch().type());

            // Count the total size we are going to write from this process
            //outputSize_ += pf1.size()*pType::nComponents*sizeof(ioScalar); // size of data array

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

        Info<< "    fieldWrite: " << varPath << endl;
        // {
        //     OBufStream os(iobuffer_, adiosCore::strFormat);
        //     os << field;
        // }

        {
            OCompactBufStream os(iobuffer_, adiosCore::strFormat);
            os << field;
        }

        // Do the actual write
        writeVariable(varPath, iobuffer_.cdata());

#ifdef FOAM_ADIOS_PATCH_WRITE

        // ADIOS requires non-const access to the field for writing (?)
        FieldType& fld = const_cast<FieldType&>(field);

        // Do the actual write
        writeVariable(datasetName, fld.cdata());
        // writeVariable(datasetName, fld.internalField().cdata());

        typename FieldType::GeometricBoundaryField& bf = fld.boundaryField();

        // Do the same for all patches
        forAll(bf, patchI)
        {
            typename FieldType::PatchFieldType& pf1 = bf[patchI];

            fileName patchName(datasetName/"patch" + Foam::name(patchI));
            Info<< "      patchfield " << patchI << ":" << endl;

#if 0
            // FIXME: what's the name of pf1 in the output?
            adios_write(fileID_, patchName.c_str(), pf1.data());
#else
            OStringStream os;
            os << pf1;
            Info<< "patch " << patchName << " content: " << os.str() << endl;
            writeVariable(patchName, os);
#endif
        }

#endif /* FOAM_ADIOS_PATCH_WRITE */

    }
}


// ************************************************************************* //
