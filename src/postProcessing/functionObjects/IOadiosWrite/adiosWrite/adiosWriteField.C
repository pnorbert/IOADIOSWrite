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
#include "basicKinematicCloud.H"
#include "nullObject.H"
#include "FlatListOutput.H"
#include "IOstreams.H"
#include "pointFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::adiosWrite::regionInfo::classifyFields
(
    const fvMesh& mesh,
    bool verbose
)
{
    if (verbose)
    {
        Info<< "  " << info() << endl;
    }

    fieldsToWrite_.clear(); // reset

    // clouds are handled elsewhere:
    wordHashSet ignore(mesh.names<cloud>());
    HashTable<word> unsupported;

    const wordList allFields = mesh.names();
    forAll(allFields, i)
    {
        const word& name = allFields[i];
        if (ignore.found(name) || findStrings(ignoredFields_, name))
        {
            continue;
        }

        const regIOobject* obj = mesh.find(name)();
        const word& type = obj->type();

        // auto-write field or explicitly requested
        if
        (
            (autoWrite() && obj->writeOpt() == IOobject::AUTO_WRITE)
         || findStrings(requestedFields_, name)
        )
        {
            if (supportedFieldType(type))
            {
                fieldsToWrite_.insert(name, type);
            }
            else
            {
                unsupported.set(name, type);
            }
        }
    }

    if (verbose)
    {
        wordList names = fieldsToWrite_.sortedToc();
        forAll(names, fieldI)
        {
            Info<< "    " << names[fieldI] << "  ("
                << fieldsToWrite_[names[fieldI]] << ")" << endl;
        }
    }

    if (!unsupported.empty())
    {
        wordList names = unsupported.sortedToc();
        wordList types(names.size());

        forAll(names, fieldI)
        {
            types[fieldI] = unsupported[names[fieldI]];
        }

        WarningInFunction
            << nl
            << unsupported.size() << " fields not handled by adiosWrite" << nl
            << "  names: " << FlatListOutput<word>(names) << nl
            << "  types: " << FlatListOutput<word>(types) << nl << endl;
    }


    return fieldsToWrite_.size();
}


Foam::label Foam::adiosWrite::fieldDefine(const regionInfo& rInfo)
{
    label nFields = 0;
    Info<< "  adiosWrite::fieldDefine: " << rInfo.info() << endl;

    const fvMesh&  mesh = time_.lookupObject<fvMesh>(rInfo.name());
    const fileName path = rInfo.fieldPath();
    wordList fieldNames = rInfo.fieldsToWrite_.sortedToc();

    forAll(fieldNames, fieldI)
    {
        const word& fieldName = fieldNames[fieldI];
        const word& fieldType = rInfo.fieldsToWrite_[fieldName];

        if (!(static_cast<const objectRegistry&>(mesh).found(fieldName)))
        {
            continue;
        }

        regIOobject* obj = mesh.find(fieldName)();

        int64_t varid = -1;

        // point fields
        if (fieldType == pointScalarField::typeName)
        {
            varid = defineField
            (
                static_cast<pointScalarField&>(*obj),
                path
            );
        }
        else if (fieldType == pointVectorField::typeName)
        {
            varid = defineField
            (
                static_cast<pointVectorField&>(*obj),
                path
            );
        }
        else if (fieldType == pointSphericalTensorField::typeName)
        {
            varid = defineField
            (
                static_cast<pointSphericalTensorField&>(*obj),
                path
            );
        }
        else if (fieldType == pointSymmTensorField::typeName)
        {
            varid = defineField
            (
                static_cast<pointSymmTensorField&>(*obj),
                path
            );
        }
        else if (fieldType == pointTensorField::typeName)
        {
            varid = defineField
            (
                static_cast<pointTensorField&>(*obj),
                path
            );
        }

        // surface fields
        else if (fieldType == surfaceScalarField::typeName)
        {
            varid = defineField
            (
                static_cast<surfaceScalarField&>(*obj),
                path
            );
        }
        else if (fieldType == surfaceVectorField::typeName)
        {
            varid = defineField
            (
                static_cast<surfaceVectorField&>(*obj),
                path
            );
        }
        else if (fieldType == surfaceSphericalTensorField::typeName)
        {
            varid = defineField
            (
                static_cast<surfaceSphericalTensorField&>(*obj),
                path
            );
        }
        else if (fieldType == surfaceSymmTensorField::typeName)
        {
            varid = defineField
            (
                static_cast<surfaceSymmTensorField&>(*obj),
                path
            );
        }
        else if (fieldType == surfaceTensorField::typeName)
        {
            varid = defineField
            (
                static_cast<surfaceTensorField&>(*obj),
                path
            );
        }

        // volume fields
        else if (fieldType == volScalarField::typeName)
        {
            varid = defineField
            (
                static_cast<volScalarField&>(*obj),
                path
            );
        }
        else if (fieldType == volVectorField::typeName)
        {
            varid = defineField
            (
                static_cast<volVectorField&>(*obj),
                path
            );
        }
        else if (fieldType == volSphericalTensorField::typeName)
        {
            varid = defineField
            (
                static_cast<volSphericalTensorField&>(*obj),
                path
            );
        }
        else if (fieldType == volSymmTensorField::typeName)
        {
            varid = defineField
            (
                static_cast<volSymmTensorField&>(*obj),
                path
            );
        }
        else if (fieldType == volTensorField::typeName)
        {
            varid = defineField
            (
                static_cast<volTensorField&>(*obj),
                path
            );
        }

        // internal volume fields
        else if (fieldType == volScalarField::DimensionedInternalField::typeName)
        {
            varid = defineInternalField
            (
                static_cast<volScalarField::DimensionedInternalField&>(*obj),
                path
            );
        }
        else if (fieldType == volVectorField::DimensionedInternalField::typeName)
        {
            varid = defineInternalField
            (
                static_cast<volVectorField::DimensionedInternalField&>(*obj),
                path
            );
        }

        if (varid != -1)
        {
            ++nFields;
        }
    }

    return nFields;
}


void Foam::adiosWrite::fieldWrite(const regionInfo& rInfo)
{
    Info<< "  adiosWrite::fieldWrite: " << rInfo.info() << endl;

    const fvMesh&  mesh = time_.lookupObject<fvMesh>(rInfo.name());
    const fileName path = rInfo.fieldPath();
    wordList fieldNames = rInfo.fieldsToWrite_.sortedToc();

    forAll(fieldNames, fieldI)
    {
        const word& fieldName = fieldNames[fieldI];
        const word& fieldType = rInfo.fieldsToWrite_[fieldName];

        if (!static_cast<const objectRegistry&>(mesh).found(fieldName))
        {
            continue;
        }

        regIOobject* obj = mesh.find(fieldName)();

        // point fields
        if (fieldType == pointScalarField::typeName)
        {
            writeField
            (
                static_cast<pointScalarField&>(*obj),
                path
            );
        }
        else if (fieldType == pointVectorField::typeName)
        {
            writeField
            (
                static_cast<pointVectorField&>(*obj),
                path
            );
        }
        else if (fieldType == pointSphericalTensorField::typeName)
        {
            writeField
            (
                static_cast<pointSphericalTensorField&>(*obj),
                path
            );
        }
        else if (fieldType == pointSymmTensorField::typeName)
        {
            writeField
            (
                static_cast<pointSymmTensorField&>(*obj),
                path
            );
        }
        else if (fieldType == pointTensorField::typeName)
        {
            writeField
            (
                static_cast<pointTensorField&>(*obj),
                path
            );
        }

        // surface fields
        else if (fieldType == surfaceScalarField::typeName)
        {
            writeField
            (
                static_cast<surfaceScalarField&>(*obj),
                path
            );
        }
        else if (fieldType == surfaceVectorField::typeName)
        {
            writeField
            (
                static_cast<surfaceVectorField&>(*obj),
                path
            );
        }
        else if (fieldType == surfaceSphericalTensorField::typeName)
        {
            writeField
            (
                static_cast<surfaceSphericalTensorField&>(*obj),
                path
            );
        }
        else if (fieldType == surfaceSymmTensorField::typeName)
        {
            writeField
            (
                static_cast<surfaceSymmTensorField&>(*obj),
                path
            );
        }
        else if (fieldType == surfaceTensorField::typeName)
        {
            writeField
            (
                static_cast<surfaceTensorField&>(*obj),
                path
            );
        }

        // volume fields
        else if (fieldType == volScalarField::typeName)
        {
            writeField
            (
                static_cast<volScalarField&>(*obj),
                path
            );
        }
        else if (fieldType == volVectorField::typeName)
        {
            writeField
            (
                static_cast<volVectorField&>(*obj),
                path
            );
        }
        else if (fieldType == volSphericalTensorField::typeName)
        {
            writeField
            (
                static_cast<volSphericalTensorField&>(*obj),
                path
            );
        }
        else if (fieldType == volSymmTensorField::typeName)
        {
            writeField
            (
                static_cast<volSymmTensorField&>(*obj),
                path
            );
        }
        else if (fieldType == volTensorField::typeName)
        {
            writeField
            (
                static_cast<volTensorField&>(*obj),
                path
            );
        }

        // internal volume fields
        else if (fieldType == volScalarField::DimensionedInternalField::typeName)
        {
            writeField
            (
                static_cast<volScalarField::DimensionedInternalField&>(*obj),
                path
            );
        }
        else if (fieldType == volVectorField::DimensionedInternalField::typeName)
        {
            writeField
            (
                static_cast<volVectorField::DimensionedInternalField&>(*obj),
                path
            );
        }
    }
}


// ************************************************************************* //
