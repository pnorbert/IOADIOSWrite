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
    OutputCounter os(adiosCore::strFormat);

    size_t bufLen = 0;
    size_t maxLen = 0;

    forAll(fields, fieldI)
    {
        // Lookup field
        const FieldType& field = mesh.lookupObject<FieldType>(fields[fieldI]);
        const typename FieldType::GeometricBoundaryField& bfield =
            field.boundaryField();

        fileName varPath = rInfo.fieldVarPath(fields[fieldI]);

        os.rewind();
        os << field;
        bufLen = os.size();
        maxLen = Foam::max(maxLen, bufLen);

//         OutputStringStreamer check(adiosCore::strFormat);
//         check << field;
//         size_t oslen = check.str().size();
//
//         Pout<< "    fieldDefine: " << varPath
//             << "  (size " << field.size() << ")"
//             << "  stream-size " << oslen << " (counted " << bufLen << ")" << endl
//             << "  stream-content " << check.str().c_str() << endl
//             << "  ----" << endl;

        defineVariable(varPath, adios_unsigned_byte, bufLen);

        // volScalarField etc.
        defineAttribute("class", varPath, field.type());

        // could also write field.dimensions() as an attribute
        // if needed to save parsing

        // independent of how we store fields,
        // a quick lookup of field patch types could prove useful
        forAll(bfield, patchI)
        {
            const typename FieldType::PatchFieldType& pf = bfield[patchI];
            fileName patchPath = varPath/"patch" + Foam::name(patchI);

            defineAttribute("type", patchPath, pf.type());
        }


#ifdef FOAM_ADIOS_PATCH_WRITE
        varPath = "region" + Foam::name(rInfo.index_) / "fields" / fields[fieldI];

        string localDim = Foam::name(field.size());
        const int nCmpt = pTraits<typename FieldType::value_type>::nComponents;

        if (nCmpt > 1)
        {
            localDim += "," + Foam::name(nCmpt);
        }

        // store local to this process
        // Type is float or double depending on OpenFOAM precision
        adios_define_var
        (
            groupID_,
            varPath.c_str(),
            "",
            ADIOS_SCALAR,
            localDim.c_str(),
            "",
            ""
        );

        // Count the total size we are going to write from this process
        outputSize_ += nCmpt * field.size() * sizeof(ioScalar);

        forAll(bfield, patchI)
        {
            const typename FieldType::PatchFieldType& pf = bfield[patchI];
            fileName patchPath = varPath/"patch" + Foam::name(patchI);

            Info<< "      patchfield " << patchI
                << "(" << fields[fieldI] << ")"
                << ": name=" << pf.patch().name()
                << ": type=" << pf.type()
                << " empty=" << pf.empty()
                << " size=" << pf.size()
                << endl;

            os.rewind();
            os << pf;
            bufLen = os.size();
            maxLen = Foam::max(maxLen, bufLen);

            defineVariable(patchPath, adios_unsigned_byte, bufLen);

            // additional attributes to describe this patch
            defineAttribute("name", patchPath, pf.patch().name());
            defineAttribute("type", patchPath, pf.type());

            // Need much more such logic if we want to store patch values as arrays.
            // This is non-trivial since not every patch-field has values for it faces

            // Count the total size we are going to write from this process
            // outputSize_ += pf.size()*pType::nComponents*sizeof(ioScalar); // size of data array

            // forAll(pf, faceI)
            // {
            //     pf[faceI] = ...
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
    forAll(fields, fieldI)
    {
        // Lookup field
        const FieldType& field = mesh.lookupObject<FieldType>(fields[fieldI]);

        fileName varPath = rInfo.fieldVarPath(fields[fieldI]);

        Info<< "    fieldWrite: " << varPath << endl;
        {
            OutputBufStreamer os(iobuffer_, adiosCore::strFormat);
            os << field;
        }

        // Do the actual write (as stream)
        writeVariable(varPath, iobuffer_.cdata());

#ifdef FOAM_ADIOS_PATCH_WRITE
        varPath = "region" + Foam::name(rInfo.index_) / "fields" / fields[fieldI];

        Info<< "    fieldWrite: " << varPath << endl;

        // write internalField
        writeVariable(varPath, field.internalField().cdata());

        const typename FieldType::GeometricBoundaryField& bfield =
            field.boundaryField();

        // Do the same for all patches
        forAll(bfield, patchI)
        {
            const typename FieldType::PatchFieldType& pf = bfield[patchI];

            fileName patchPath = varPath / "patch" + Foam::name(patchI);
            // Info<< "      patchfield " << patchI << ":" << endl;

            {
                OutputBufStreamer os(iobuffer_, adiosCore::strFormat);
                os << pf;
            }

            writeVariable(patchPath, iobuffer_.cdata());
        }

#endif /* FOAM_ADIOS_PATCH_WRITE */

    }
}


// ************************************************************************* //
