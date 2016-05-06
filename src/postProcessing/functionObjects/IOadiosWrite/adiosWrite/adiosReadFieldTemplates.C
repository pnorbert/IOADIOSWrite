/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
     \\/     M anipulation  |
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
#include "adiosReader.H"

#include "dictionary.H"

#include "IBufStream.H"
#include "IOstream.H"
#include "Ostream.H"
#include "OFstream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FieldType>
bool Foam::adiosWrite::fieldRead
(
    FieldType& field,
    adiosReader::helper& helper,
    const adiosReader::fieldInfo& src
)
{
    bool ok = true;

    // Read data from file - fatal error if this fails!
    ok = helper.getDataSet(src);
    if (ok)
    {
        // read fields via dictionary
        IBufStream is(helper.buffer, adiosCore::strFormat);
        dictionary dict(is);

        Pout<<"dictionary: " << field.name() << " with "
            << dict.toc() << " boundaryField: "
            << dict.subDict("boundaryField").toc() << endl;

        // could also verify dimensions
        field.readField(dict, "internalField");
        field.boundaryField().readField(field, dict.subDict("boundaryField"));

        // TODO: adjust for referenceLevel?
        //
        //     if (dict.found("referenceLevel"))
        //     {
        //         Type fieldAverage(pTraits<Type>(dict.lookup("referenceLevel")));
        //         Field<Type>::operator+=(fieldAverage);
        //         forAll(boundaryField_, patchI)
        //         {
        //             boundaryField_[patchI] == boundaryField_[patchI] + fieldAverage;
        //         }
        //     }
    }

    return ok;
}


template<class FieldType>
bool Foam::adiosWrite::fieldRead
(
    adiosReader::helper& helper,
    const fvMesh& mesh,
    const adiosReader::fieldInfo& src
)
{
    bool ok = true;

    // Lookup field
    FieldType& field = const_cast<FieldType&>
    (
        mesh.lookupObject<FieldType>(src.name())
    );

    Pout<< "    readField via dictionary: " << field.name() << endl;

    // Read data from file - fatal error if this fails!
    ok = helper.getDataSet(src);
    if (ok)
    {
        // read fields via dictionary
        IBufStream is(helper.buffer, adiosCore::strFormat);
        dictionary dict(is);

        Pout<<"dictionary: " << field.name() << " with "
            << dict.toc() << " boundaryField: "
            << dict.subDict("boundaryField").toc() << endl;

        // could also verify dimensions
        field.readField(dict, "internalField");
        field.boundaryField().readField(field, dict.subDict("boundaryField"));

        // TODO: adjust for referenceLevel?
        //
        //     if (dict.found("referenceLevel"))
        //     {
        //         Type fieldAverage(pTraits<Type>(dict.lookup("referenceLevel")));
        //         Field<Type>::operator+=(fieldAverage);
        //         forAll(boundaryField_, patchI)
        //         {
        //             boundaryField_[patchI] == boundaryField_[patchI] + fieldAverage;
        //         }
        //     }
    }

    return ok;
}


// ************************************************************************* //
