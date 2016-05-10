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

#include "ParticleBinaryBlob.H"

#include "particle.H"
#include "IStringStream.H"
#include "token.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ParticleBinaryBlob::clear()
{
    static_cast<Container&>(*this).clear();
    types_.clear();
    names_.clear();
}


void Foam::ParticleBinaryBlob::setTypes
(
    const UList<word>& inputTypes,
    const bool raw
)
{
    clear();
    size_t offset = (raw ? 0 : 1);

    forAll(inputTypes, elemI)
    {
        Fragment frag(inputTypes[elemI], offset);
        offset = frag.end();

        // This is a hack, but we are presented with at least two blocks
        // of binary content, each of which is surrounded by '()'
        // content1 => (particle)
        // content2 => (other)

        if (!raw && offset == Foam::particle::sizeofFields + 1)
        {
            // 2 == closing ')' from particle + opening '(' for next sequence
            offset += 2;
        }

        append(frag);
    }
}


void Foam::ParticleBinaryBlob::setTypes
(
    Istream& is,
    const bool raw
)
{
    DynamicList<word> lst(32);

    token tok;

    // the input types are the driver
    while (is.good())
    {
        is.read(tok);

        if (tok.type() != token::WORD)
        {
            // or warn and ignore the balance of the input
            FatalErrorInFunction
                << "Unsupported particle PropertyType '"
                << tok.type() << "'"
                << exit(FatalError);
        }
        else
        {
            lst.append(tok.wordToken());
        }
    }

    setTypes(lst, raw);
}


void Foam::ParticleBinaryBlob::setTypes
(
    const string& inputTypes,
    const bool raw
)
{
    IStringStream is(inputTypes);
    setTypes(is, raw);
}


void Foam::ParticleBinaryBlob::setNames
(
    const UList<word>& inputNames
)
{
    names_.clear();

    if (inputNames.size() != size())
    {
        FatalErrorInFunction
            << "Mismatch in number of names ("
            << inputNames.size()
            << ") and number of particle-binary-blob fragments ("
            << size() << ")"
            << exit(FatalError);
    }

    // even more safety
    label elemI = 0;
    const label nNames = Foam::min(inputNames.size(), size());

    forAllIter(Container, *this, iter)
    {
        if (elemI < nNames)
        {
            iter().name_ = inputNames[elemI++];
        }
        else
        {
            break;
        }
    }
}


void Foam::ParticleBinaryBlob::setNames
(
    Istream& is
)
{
    names_.clear();

    token tok;

    // the input types are the driver, get the corresponding name
    label nNames = 0;
    forAllIter(Container, *this, iter)
    {
        Fragment& frag = iter();
        // Info<<"filling: " << frag.type_ << endl;

        is.read(tok);

        if (tok.type() == token::WORD)
        {
            // name -> use directly
            frag.name_ = tok.wordToken();
            ++nNames;
        }
        else if
        (
            tok.type() == token::PUNCTUATION
         && tok.pToken() == token::BEGIN_LIST
        )
        {
            // should be a list of names.
            // Eg, "(UTurbx UTurby UTurbz)"

            DynamicList<word> cmptNames;

            bool readingList = true;
            while (readingList && is.good())
            {
                is.read(tok);

                if (tok.type() == token::WORD)
                {
                    cmptNames.append(tok.wordToken());
                }
                else if
                (
                    tok.type() == token::PUNCTUATION
                 && tok.pToken() == token::END_LIST
                )
                {
                    readingList = false;
                }
            }

            // reduce list of names:
            word varName;
            if (cmptNames.size() == 3 && frag.nComponents == 3)
            {
                // reduce things like "(Ux Uy Uz)" -> "U"
                // builtin aliases:
                // - "(Px Py Pz)" -> "position"

                const word& cmpt0 = cmptNames[0];
                if (cmpt0 == "Px")
                {
                    varName = "position";
                }
                else if (cmpt0[cmpt0.size()-1] == 'x')
                {
                    varName = cmpt0.substr(0, cmpt0.size()-1);
                }
            }

            if (varName.empty())
            {
                // previous heuristics failed,
                // simply concatenate with '-' delimiter
                forAll(cmptNames, elemI)
                {
                    if (elemI > 0)
                    {
                        varName += '-';
                    }
                    varName += cmptNames[elemI];
                }
            }

            if (varName.empty())
            {
                // TODO: terminate here with the rest being treated
                // as trailing bytes.
                // For now, just use leave as 'unnamed'
                break;
            }
            else
            {
                frag.name_ = varName;
            }

            ++nNames;
        }
        else
        {
            // FATAL!
            Info<< "bad-input: " << tok << endl;
        }
    }

    // Info<<"read in " << size() << " types and " << nNames << " names" << endl;
}


void Foam::ParticleBinaryBlob::setNames
(
    const string& inputNames
)
{
    IStringStream is(inputNames);
    setNames(is);
}


void Foam::ParticleBinaryBlob::makeSummary()
{
    types_.clear();
    names_.clear();

    types_.setSize(size());
    names_.setSize(size());

    label idx = 0;
    forAllConstIter(Container, *this, iter)
    {
        types_[idx] = iter().type_;
        names_[idx] = iter().name_;

        ++idx;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ParticleBinaryBlob::ParticleBinaryBlob()
:
    Container(),
    types_(),
    names_()
{}



Foam::ParticleBinaryBlob::ParticleBinaryBlob
(
    const UList<word>& inputTypes,
    const UList<word>& inputNames,
    const UList<int>& offsetList,
    const UList<int>& bytesList
)
:
    Container(),
    types_(),
    names_()
{
    const label size = inputTypes.size();

    if
    (
        size != inputNames.size()
     || size != offsetList.size()
     || size != bytesList.size()
    )
    {
        FatalErrorInFunction
            << "Mismatch in particle-binary-blob sizes (types:"
            << size
            << " names:" << inputTypes.size()
            << " offsets:" << offsetList.size()
            << " bytes:" << bytesList.size() << ")"
            << exit(FatalError);
    }

    forAll(inputTypes, elemI)
    {
        Fragment frag(inputTypes[elemI], offsetList[elemI]);
        frag.name_  = inputNames[elemI];
        frag.width_ = bytesList[elemI];

        append(frag);
    }

    makeSummary();
}


Foam::ParticleBinaryBlob::ParticleBinaryBlob
(
    const UList<word>& inputTypes,
    const UList<word>& inputNames,
    const bool raw
)
:
    Container(),
    types_(),
    names_()
{
    setTypes(inputTypes, raw);
    setNames(inputNames);

    makeSummary();
}


Foam::ParticleBinaryBlob::ParticleBinaryBlob
(
    const string& inputTypes,
    const string& inputNames,
    const bool raw
)
:
    Container(),
    types_(),
    names_()
{
    setTypes(inputTypes, raw);
    setNames(inputNames);

    makeSummary();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::wordList& Foam::ParticleBinaryBlob::types() const
{
    return types_;
}


const Foam::wordList& Foam::ParticleBinaryBlob::names() const
{
    return names_;
}


Foam::List<int> Foam::ParticleBinaryBlob::offsets() const
{
    List<int> lst(size());

    label i=0;
    forAllConstIter(Container, *this, iter)
    {
        lst[i++] = iter().offset();
    }

    return lst;
}


Foam::List<int> Foam::ParticleBinaryBlob::sizes() const
{
    List<int> lst(size());

    label i=0;
    forAllConstIter(Container, *this, iter)
    {
        lst[i++] = iter().sizeOf();
    }

    return lst;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ParticleBinaryBlob& blob
)
{
    // Write contents
    os << static_cast<const ParticleBinaryBlob::Container&>(blob) << nl;

    os.check("Ostream& operator<<(Ostream&, const ParticleBinaryBlob&)");

    return os;
}


// ************************************************************************* //
