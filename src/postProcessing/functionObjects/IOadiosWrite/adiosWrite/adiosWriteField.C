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
    if (autoWrite())
    {
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
                obj->writeOpt() == IOobject::AUTO_WRITE
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
    }
    else
    {
        labelList indices = findStrings(requestedFields_, allFields);

        forAll(indices, fieldI)
        {
            const word& name = allFields[indices[fieldI]];
            if (ignore.found(name) || findStrings(ignoredFields_, name))
            {
                continue;
            }

            const regIOobject* obj = mesh.find(name)();
            const word& type = obj->type();

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


Foam::label Foam::adiosWrite::classifyFields(bool verbose)
{
    if (verbose)
    {
        Info<< endl << "Foam::adiosWrite::classifyFields:" << endl;
    }

    label nFields = 0;
    forAllIter(RegionInfoContainer, regions_, iter)
    {
        regionInfo& r = iter();

        nFields += r.classifyFields
        (
            time_.lookupObject<fvMesh>(r.name()),
            verbose
        );
    }

    return nFields;
}


Foam::label Foam::adiosWrite::fieldDefine(const regionInfo& rInfo)
{
    label nFields = 0;
    Info<< "  adiosWrite::fieldDefine: " << rInfo.info() << endl;

    const fvMesh& mesh = time_.lookupObject<fvMesh>(rInfo.name());
    wordList fieldNames = rInfo.fieldsToWrite_.sortedToc();
    forAll(fieldNames, fieldI)
    {
        const word& fieldName = fieldNames[fieldI];
        const word& fieldType = rInfo.fieldsToWrite_[fieldName];
        const fileName varPath = rInfo.fieldPath(fieldName);

        if (!(static_cast<const objectRegistry&>(mesh).found(fieldName)))
        {
            continue;
        }

        regIOobject* obj = mesh.find(fieldName)();

        int64_t varid = -1;

        if (fieldType == volScalarField::typeName)
        {
            varid = defineField
            (
                static_cast<volScalarField&>(*obj),
                varPath
            );
        }
        else if (fieldType == volVectorField::typeName)
        {
            varid = defineField
            (
                static_cast<volVectorField&>(*obj),
                varPath
            );
        }
        else if (fieldType == surfaceScalarField::typeName)
        {
            varid = defineField
            (
                static_cast<surfaceScalarField&>(*obj),
                varPath
            );
        }
        else if (fieldType == volSphericalTensorField::typeName)
        {
            varid = defineField
            (
                static_cast<volSphericalTensorField&>(*obj),
                varPath
            );
        }
        else if (fieldType == volSymmTensorField::typeName)
        {
            varid = defineField
            (
                static_cast<volSymmTensorField&>(*obj),
                varPath
            );
        }
        else if (fieldType == volTensorField::typeName)
        {
            varid = defineField
            (
                static_cast<volTensorField&>(*obj),
                varPath
            );
        }
        else if (fieldType == volScalarField::DimensionedInternalField::typeName)
        {
            varid = defineInternalField
            (
                static_cast<volScalarField::DimensionedInternalField&>(*obj),
                varPath
            );
        }
        else if (fieldType == volVectorField::DimensionedInternalField::typeName)
        {
            varid = defineInternalField
            (
                static_cast<volVectorField::DimensionedInternalField&>(*obj),
                varPath
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

    const fvMesh& mesh = time_.lookupObject<fvMesh>(rInfo.name());

    wordList fieldNames = rInfo.fieldsToWrite_.sortedToc();
    forAll(fieldNames, fieldI)
    {
        const word& fieldName = fieldNames[fieldI];
        const word& fieldType = rInfo.fieldsToWrite_[fieldName];
        const fileName varPath = rInfo.fieldPath(fieldName);

        if (!static_cast<const objectRegistry&>(mesh).found(fieldName))
        {
            continue;
        }

        regIOobject* obj = mesh.find(fieldName)();

        if (fieldType == volScalarField::typeName)
        {
            writeField
            (
                static_cast<volScalarField&>(*obj),
                varPath
            );
        }
        else if (fieldType == volVectorField::typeName)
        {
            writeField
            (
                static_cast<volVectorField&>(*obj),
                varPath
            );
        }
        else if (fieldType == surfaceScalarField::typeName)
        {
            writeField
            (
                static_cast<surfaceScalarField&>(*obj),
                varPath
            );
        }
        else if (fieldType == volSphericalTensorField::typeName)
        {
            writeField
            (
                static_cast<volSphericalTensorField&>(*obj),
                varPath
            );
        }
        else if (fieldType == volSymmTensorField::typeName)
        {
            writeField
            (
                static_cast<volSymmTensorField&>(*obj),
                varPath
            );
        }
        else if (fieldType == volTensorField::typeName)
        {
            writeField
            (
                static_cast<volTensorField&>(*obj),
                varPath
            );
        }
        else if (fieldType == volScalarField::DimensionedInternalField::typeName)
        {
            writeField
            (
                static_cast<volScalarField::DimensionedInternalField&>(*obj),
                varPath
            );
        }
        else if (fieldType == volVectorField::DimensionedInternalField::typeName)
        {
            writeField
            (
                static_cast<volVectorField::DimensionedInternalField&>(*obj),
                varPath
            );
        }
    }
}


// ************************************************************************* //
