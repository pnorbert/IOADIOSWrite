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
#include "nullObject.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


bool Foam::adiosWrite::supportedFieldType
(
    const word& fieldType
)
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
     // || fieldType == volSphericalTensorField::typeName
     // || fieldType == volSymmTensorField::typeName
     // || fieldType == volTensorField::typeName
    );
}


Foam::label Foam::adiosWrite::regionInfo::appendFieldGroup
(
    const word& fieldName,
    const word& fieldType
)
{
    if (fieldType == volScalarField::typeName)
    {
        scalarFields_.append(fieldName);
        return 1;
    }
    else if (fieldType == volVectorField::typeName)
    {
        vectorFields_.append(fieldName);
        return 1;
    }
    else if (fieldType == surfaceScalarField::typeName)
    {
        surfaceScalarFields_.append(fieldName);
        return 1;
    }
    /*
    else if (fieldType == volSphericalTensorField::typeName)
    {
        sphericalTensorFields_.append(fieldName);
        return 1;
    }
    else if (fieldType == volSymmTensorField::typeName)
    {
        symmTensorFields_.append(fieldName);
        return 1;
    }
    else if (fieldType == volTensorField::typeName)
    {
        tensorFields_.append(fieldName);
        return 1;
    }
    */
    else
    {
        WarningInFunction
            << "Field type " << fieldType
            << " of the field " << fieldName
            << " is not handled by adiosWrite."
            << endl;
    }

    return 0;
}


Foam::label Foam::adiosWrite::regionInfo::classifyFields
(
    const fvMesh& mesh
)
{
    label nFields = 0;
    wordList allFields = mesh.sortedNames();

    Info<< "  " << info() << endl;

    clearFields(); // clear it because we will add all of them again and again

    if (autoWrite())
    {
        forAll(allFields, i)
        {
            const word& name = allFields[i];
            const regIOobject* obj = mesh.find(name)();
            const word& type = obj->type();

            // auto-write field or explicitly requested
            bool shouldWrite =
            (
                obj->writeOpt() == IOobject::AUTO_WRITE
             || findStrings(objectNames_, name)
            );

            if (shouldWrite)
            {
                Info<< "    name = " << name << " type = " << type << endl;
                nFields += appendFieldGroup(name, type);
            }
        }
    }
    else
    {
        labelList indices = findStrings(objectNames_, allFields);

        forAll(indices, fieldI)
        {
            const word& name = allFields[indices[fieldI]];
            const word& type = mesh.find(name)()->type();

            Info<< "    name = " << name << " type = " << type << endl;
            nFields += appendFieldGroup(name, type);
        }
    }

    return nFields;
}


Foam::label Foam::adiosWrite::classifyFields()
{
    Info<< endl << "Foam::adiosWrite::classifyFields: " << endl;

    label nFields = 0;
    forAll(regions_, i)
    {
        regionInfo& r = regions_[i];

        nFields += r.classifyFields(time_.lookupObject<fvMesh>(r.name_));
    }

    return nFields;
}


// ************************************************************************* //
