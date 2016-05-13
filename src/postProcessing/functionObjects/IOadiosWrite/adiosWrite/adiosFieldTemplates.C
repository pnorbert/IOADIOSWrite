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
#include "IOstreams.H"
#include "Ostream.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


template<class FieldType>
bool Foam::adiosWrite::fieldRead
(
    FieldType& field,
    const adiosReader& reader,
    const adiosReader::fieldInfo& src
)
{
    // Read data from file - fatal error if this fails!
    size_t nread = reader.getBuffered(src.fullName(), iobuffer_);
    if (nread)
    {
        // read fields via dictionary
        IBufStream is(iobuffer_, nread, adiosCore::strFormat);
        dictionary dict(is);

        // Pout<<"dictionary: " << field.name() << " with "
        //     << dict.toc() << " boundaryField: "
        //     << dict.subDict("boundaryField").toc() << endl;

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

    return nread;
}


template<class FieldType>
bool Foam::adiosWrite::fieldRead
(
    const adiosReader& reader,
    const fvMesh& mesh,
    const adiosReader::fieldInfo& src
)
{
    // Lookup field
    FieldType& field = const_cast<FieldType&>
    (
        mesh.lookupObject<FieldType>(src.name())
    );

    // Pout<< "    readField via dictionary: " << field.name() << endl;

    // Read data from file - fatal error if this fails!
    size_t nread = reader.getBuffered(src.fullName(), iobuffer_);

    if (nread)
    {
        // read fields via dictionary
        IBufStream is(iobuffer_, nread, adiosCore::strFormat);
        dictionary dict(is);

        // Pout<<"dictionary: " << field.name() << " with "
        //     << dict.toc() << " boundaryField: "
        //     << dict.subDict("boundaryField").toc() << endl;

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

    return nread;
}


// ************************************************************************* //
