/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 Norbert Podhorszki
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

        const fileName varPath = rInfo.fieldPath(fields[fieldI]);

        // volScalarField etc.
        defineAttribute("class", varPath, field.type());

        // could also write field.dimensions() as an attribute
        // if needed to save parsing

        // independent of how we store fields,
        // a quick lookup of field patch types may prove useful
        stringList pTypes(bfield.size());

        forAll(bfield, patchI)
        {
            const typename FieldType::PatchFieldType& pf = bfield[patchI];

            pTypes[patchI] = pf.type();
        }

        defineListAttribute("patch-types", varPath, pTypes);

        os.rewind();
        os << field;

        bufLen = defineStreamVariable(varPath, os.size());
        maxLen = Foam::max(maxLen, bufLen);

    }

    // Pout<< "max stream-size: " << maxLen << endl;

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

        const fileName varPath = rInfo.fieldPath(fields[fieldI]);
        Info<< "    fieldWrite: " << varPath << endl;
        {
            OutputBufStreamer os(iobuffer_, adiosCore::strFormat);
            os << field;
        }

        // Do the actual write (as stream)
        writeVariable(varPath, iobuffer_);
    }
}


// ************************************************************************* //
