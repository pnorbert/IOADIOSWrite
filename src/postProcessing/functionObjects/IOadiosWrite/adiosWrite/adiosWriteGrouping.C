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


Foam::label Foam::adiosWrite::appendFieldGroup
(
    regionInfo& r,
    const word& fieldName,
    const word& fieldType
)
{
    if (fieldType == volScalarField::typeName)
    {
        r.scalarFields_.append(fieldName);
        return 1;
    }
    else if (fieldType == volVectorField::typeName)
    {
        r.vectorFields_.append(fieldName);
        return 1;
    }
    else if (fieldType == surfaceScalarField::typeName)
    {
        r.surfaceScalarFields_.append(fieldName);
        return 1;
    }
    /*
    else if (fieldType == volSphericalTensorField::typeName)
    {
        //r.sphericalTensorFields_.append(fieldName);
        return 0;
    }
    else if (fieldType == volSymmTensorField::typeName)
    {
        //r.symmTensorFields_.append(fieldName);
        return 0;
    }
    else if (fieldType == volTensorField::typeName)
    {
        //r.tensorFields_.append(fieldName);
        return 0;
    }
    */
    else
    {
        WarningInFunction
            << "Field type " << fieldType
            << "of the field " << fieldName
            << " is not handled by adiosWrite."
            << endl;
    }

    return 0;
}

/*void Foam::adiosWrite::test_print_obr()
{
    Info<< "adiosWrite objectRegistry list: " << endl;
    wordList allFields = obr_.sortedNames();
    forAll(regions_, i)
    {
        regionInfo& r = regions_[i];
        labelList indices = findStrings(r.objectNames_, allFields);
        forAll(indices, fieldI)
        {
            const word& name = allFields[indices[fieldI]];
            const word& type = obr_.find(name)()->type();
            Info<< "  name = " << name << "  type = " << type << endl;

        }
    }
}*/

Foam::label Foam::adiosWrite::classifyFields()
{
    label nFields = 0;

    // test_print_obr();

    Info<< endl << "Foam::adiosWrite::classifyFields: " << endl;
    // Check currently available fields

    forAll(regions_, regionI)
    {
        regionInfo& r = regions_[regionI];
        Info<< "  region " << regionI << " " << r.name_ << ": " << endl;

        r.scalarFields_.clear(); // clear it because we will add all of them again and again
        r.vectorFields_.clear();
        r.surfaceScalarFields_.clear();

        const fvMesh& mesh = time_.lookupObject<fvMesh>(r.name_);
        wordList allFields = mesh.sortedNames();
        labelList indices = findStrings(r.objectNames_, allFields);

        forAll(indices, fieldI)
        {
            const word& fieldName = allFields[indices[fieldI]];
            const word& type = mesh.find(fieldName)()->type();

            Info<< "    name = " << fieldName << "  type = " << type << endl;
            nFields += appendFieldGroup(r, fieldName, type);
        }
    }

    return nFields;
}


// ************************************************************************* //
