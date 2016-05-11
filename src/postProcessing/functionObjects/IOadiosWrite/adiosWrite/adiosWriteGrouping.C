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
#include "FlatListOutput.H"
#include "basicKinematicCloud.H"


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
            Info<< "    name = " << fieldsToWrite_[names[fieldI]]
                << " type = " << fieldsToWrite_[names[fieldI]]
                << endl;
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
    forAllIter(SLList<regionInfo>, regions_, iter)
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


// ************************************************************************* //
