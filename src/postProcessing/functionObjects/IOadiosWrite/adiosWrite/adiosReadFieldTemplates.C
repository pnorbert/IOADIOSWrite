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
#include "IStringStream.H"
#include "IStringStreamBuf.H"

#include "IOstream.H"
#include "Ostream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FieldType>
bool Foam::adiosWrite::fieldRead
(
    IStringStreamBuf& is,
    adiosReader::helper& helper,
    const fvMesh& mesh,
    const fieldGroup<typename FieldType::value_type>& fields,
    const label regionId
)
{
    bool ok = true;

    forAll(fields, fieldI)
    {
        // Lookup field
        FieldType& field = const_cast<FieldType&>
        (
            mesh.lookupObject<FieldType>(fields[fieldI])
        );

        Info<< "    readField via dictionary: " << field.name() << endl;

        fileName datasetName
        (
            "region" + Foam::name(regionId)
          / "field" / fields[fieldI] / "stream"
        );


        // Read data from file
        ok = helper.getDataSet(datasetName, is);
        if (ok)
        {
            Info<<"istream content:" << endl;
            Info<<"has " << is.stdStream().rdbuf()->in_avail() << " chars" << endl;
            is.print(Info);
            Info<< "is.good: " << is.good() << endl;
            Info<< "is.eof: " << is.eof() << endl;

            char c;
            while (is.good() && !is.eof())
            {
                is.get(c);
                Info<< char(c);
            }
            Info<< endl << "DONE" << endl;

            // read fields via dictionary
            // dictionary dict(is);

            // Info<<"dictionary: " << dict << endl;

            // field.readField(dict, "internalField");
            // field.boundaryField().readField(field, dict.subDict("boundaryField"));

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
        else
        {
            break;
        }
    }

    return ok;
}


// ************************************************************************* //
