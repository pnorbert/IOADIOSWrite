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
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::adiosWrite::supportedFieldType(const word& fieldType)
{
    if (isNull(fieldType) || fieldType.empty())
    {
        return false;
    }

    return
    (
        fieldType == volScalarField::typeName
     || fieldType == volVectorField::typeName
     || fieldType == surfaceScalarField::typeName
     || fieldType == volSphericalTensorField::typeName
     || fieldType == volSymmTensorField::typeName
     || fieldType == volTensorField::typeName

        // internal fields
     || fieldType == volScalarField::DimensionedInternalField::typeName
     || fieldType == volVectorField::DimensionedInternalField::typeName
    );
}


size_t Foam::adiosWrite::fieldDefine(const regionInfo& rInfo)
{
    Info<< "  adiosWrite::fieldDefine: " << rInfo.info() << endl;

    size_t bufLen = 0;
    size_t maxLen = 0;

    const fvMesh& mesh = time_.lookupObject<fvMesh>(rInfo.name_);

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

        if (fieldType == volScalarField::typeName)
        {
            bufLen = fieldDefine
            (
                static_cast<volScalarField&>(*obj),
                varPath
            );
        }
        else if (fieldType == volVectorField::typeName)
        {
            bufLen = fieldDefine
            (
                static_cast<volVectorField&>(*obj),
                varPath
            );
        }
        else if (fieldType == surfaceScalarField::typeName)
        {
            bufLen = fieldDefine
            (
                static_cast<surfaceScalarField&>(*obj),
                varPath
            );
        }
        else if (fieldType == volSphericalTensorField::typeName)
        {
            bufLen = fieldDefine
            (
                static_cast<volSphericalTensorField&>(*obj),
                varPath
            );
        }
        else if (fieldType == volSymmTensorField::typeName)
        {
            bufLen = fieldDefine
            (
                static_cast<volSymmTensorField&>(*obj),
                varPath
            );
        }
        else if (fieldType == volTensorField::typeName)
        {
            bufLen = fieldDefine
            (
                static_cast<volTensorField&>(*obj),
                varPath
            );
        }
        else if (fieldType == volScalarField::DimensionedInternalField::typeName)
        {
            // internal fields
            bufLen = fieldDefineInternal
            (
                static_cast<volScalarField::DimensionedInternalField&>(*obj),
                varPath
            );
        }
        else if (fieldType == volVectorField::DimensionedInternalField::typeName)
        {
            // internal fields
            bufLen = fieldDefineInternal
            (
                static_cast<volVectorField::DimensionedInternalField&>(*obj),
                varPath
            );
        }

        maxLen = Foam::max(maxLen, bufLen);
    }

    return maxLen;
}


void Foam::adiosWrite::fieldWrite(const regionInfo& rInfo)
{
    Info<< "  adiosWrite::fieldWrite: " << rInfo.info() << endl;

    const fvMesh& mesh = time_.lookupObject<fvMesh>(rInfo.name_);

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
            fieldWrite
            (
                static_cast<volScalarField&>(*obj),
                varPath
            );
        }
        else if (fieldType == volVectorField::typeName)
        {
            fieldWrite
            (
                static_cast<volVectorField&>(*obj),
                varPath
            );
        }
        else if (fieldType == surfaceScalarField::typeName)
        {
            fieldWrite
            (
                static_cast<surfaceScalarField&>(*obj),
                varPath
            );
        }
        else if (fieldType == volSphericalTensorField::typeName)
        {
            fieldWrite
            (
                static_cast<volSphericalTensorField&>(*obj),
                varPath
            );
        }
        else if (fieldType == volSymmTensorField::typeName)
        {
            fieldWrite
            (
                static_cast<volSymmTensorField&>(*obj),
                varPath
            );
        }
        else if (fieldType == volTensorField::typeName)
        {
            fieldWrite
            (
                static_cast<volTensorField&>(*obj),
                varPath
            );
        }
        else if (fieldType == volScalarField::DimensionedInternalField::typeName)
        {
            // internal fields
            fieldWrite
            (
                static_cast<volScalarField::DimensionedInternalField&>(*obj),
                varPath
            );
        }
        else if (fieldType == volVectorField::DimensionedInternalField::typeName)
        {
            // internal fields
            fieldWrite
            (
                static_cast<volVectorField::DimensionedInternalField&>(*obj),
                varPath
            );
        }
    }
}


// ************************************************************************* //
