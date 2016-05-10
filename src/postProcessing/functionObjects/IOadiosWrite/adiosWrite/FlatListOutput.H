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

Class
    Foam::FlatListOutput

Description
    Simple output adapter for list output on a single line

\*---------------------------------------------------------------------------*/

#ifndef FlatListOutput_H
#define FlatListOutput_H

#include "UList.H"
#include "token.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class FlatListOutput Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class FlatListOutput
{
    const UList<T>& list_;

    const char separator_;

public:

    // Constructors

        //- Construct with default separator
        FlatListOutput(const UList<T>& input)
        :
            list_(input),
            separator_(token::SPACE)
        {}


        //- Construct
        FlatListOutput(const char separator, const UList<T>& input)
        :
            list_(input),
            separator_(separator)
        {}


    // IOstream operators

        friend Ostream& operator<<
        (
            Ostream& os,
            const FlatListOutput<T>& flat
        )
        {
            const UList<T>& list = flat.list_;
            const char sep = flat.separator_;

            if (os.format() == IOstream::ASCII || !contiguous<T>())
            {
                bool uniform = list.size() > 1 && contiguous<T>();

                if (uniform)
                {
                    forAll(list, i)
                    {
                        if (list[i] != list[0])
                        {
                            uniform = false;
                            break;
                        }
                    }
                }

                os << list.size();
                if (uniform)
                {
                    os << token::BEGIN_BLOCK << list[0] << token::END_BLOCK;
                }
                else
                {
                    os << token::BEGIN_LIST;
                    forAll(list, i)
                    {
                        if (i) os << sep;
                        os << list[i];
                    }
                    os << token::END_LIST;
                }
            }
            else
            {
                // use standard output for binary
                os << list;
            }

            return os;
        }


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
