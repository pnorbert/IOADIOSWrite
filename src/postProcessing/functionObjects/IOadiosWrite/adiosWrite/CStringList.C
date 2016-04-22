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

#include "CStringList.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CStringList::CStringList()
:
    argc_(0),
    len_(0),
    argv_(0),
    data_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::CStringList::~CStringList()
{
    clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::CStringList::clear()
{
    argc_ = 0;
    len_  = 0;

    if (data_)
    {
        delete[] data_;
        data_ = 0;
    }
    if (argv_)
    {
        delete[] argv_;
        argv_ = 0;
    }
}


int Foam::CStringList::argc() const
{
    return argc_;
}


int Foam::CStringList::size() const
{
    return argc_;
}


size_t Foam::CStringList::length() const
{
    return len_;
}


char** Foam::CStringList::argv() const
{
    return argv_;
}


char* Foam::CStringList::data() const
{
    return data_;
}


// ************************************************************************* //
