/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "temporalAveragingItem.H"
#include "IOstreams.H"
#include "dictionaryEntry.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::temporalAveragingItem::temporalAveragingItem(Istream& is)
:
    fieldName_("unknown"),
    mean_(0),
    meanFieldName_("unknown"),
    prime2Mean_(0),
    prime2MeanFieldName_("unknown"),
    primeUPrimeMean_(0),
    primeUPrimeMeanFieldName_("unknown"),
    base_(baseType::iter),
    window_(-1.0)
{
    is.check
    (
        "Foam::functionObjects::temporalAveragingItem::temporalAveragingItem"
        "(Foam::Istream&)"
    );

    const dictionaryEntry entry(dictionary::null, is);

    fieldName_ = entry.keyword();
    entry.lookup("mean") >> mean_;
    entry.lookup("prime2Mean") >> prime2Mean_;
    entry.lookup("primeUPrimeMean") >> primeUPrimeMean_;
    base_ = baseTypeNames_[entry.lookup("base")];
    window_ = entry.lookupOrDefault<scalar>("window", -1.0);
    windowName_ = entry.lookupOrDefault<word>("windowName", "");

    meanFieldName_ = fieldName_ + EXT_MEAN;
    prime2MeanFieldName_ = fieldName_ + EXT_PRIME2MEAN;
    primeUPrimeMeanFieldName_ = fieldName_ + EXT_PRIMEUPRIMEMEAN;
    if ((window_ > 0) && (windowName_ != ""))
    {
        meanFieldName_ = meanFieldName_ + "_" + windowName_;
        prime2MeanFieldName_ = prime2MeanFieldName_ + "_" + windowName_;
        primeUPrimeMeanFieldName_ = primeUPrimeMeanFieldName_ + "_" + windowName_;
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::functionObjects::operator>>
(
    Istream& is,
    temporalAveragingItem& faItem
)
{
    is.check
    (
        "Foam::Istream& Foam::operator>>"
        "(Foam::Istream&, Foam::functionObjects::temporalAveragingItem&)"
    );

    const dictionaryEntry entry(dictionary::null, is);

    faItem.fieldName_ = entry.keyword();
    entry.lookup("mean") >> faItem.mean_;
    entry.lookup("prime2Mean") >> faItem.prime2Mean_;
    entry.lookup("primeUPrimeMean") >> faItem.primeUPrimeMean_;
    faItem.base_ = faItem.baseTypeNames_[entry.lookup("base")];
    faItem.window_ = entry.lookupOrDefault<scalar>("window", -1.0);
    faItem.windowName_ = entry.lookupOrDefault<word>("windowName", "");

    faItem.meanFieldName_ = faItem.fieldName_ + temporalAveragingItem::EXT_MEAN;
    faItem.prime2MeanFieldName_ =
        faItem.fieldName_ + temporalAveragingItem::EXT_PRIME2MEAN;
    faItem.primeUPrimeMeanFieldName_ =
        faItem.fieldName_ + temporalAveragingItem::EXT_PRIMEUPRIMEMEAN;

    if ((faItem.window_ > 0) && (faItem.windowName_ != ""))
    {
        faItem.meanFieldName_ =
            faItem.meanFieldName_ + "_" + faItem.windowName_;

        faItem.prime2MeanFieldName_ =
            faItem.prime2MeanFieldName_ + "_" + faItem.windowName_;

        faItem.primeUPrimeMeanFieldName_ =
            faItem.primeUPrimeMeanFieldName_ + "_" + faItem.windowName_;
    }
    return is;
}


Foam::Ostream& Foam::functionObjects::operator<<
(
    Ostream& os,
    const temporalAveragingItem& faItem
)
{
    os.check
    (
        "Foam::Ostream& Foam::operator<<"
        "(Foam::Ostream&, const Foam::functionObjects::temporalAveragingItem&)"
    );

    os  << faItem.fieldName_ << nl << token::BEGIN_BLOCK << nl;

    os.writeKeyword("mean")
        << faItem.mean_ << token::END_STATEMENT << nl;

    os.writeKeyword("prime2Mean")
        << faItem.prime2Mean_ << token::END_STATEMENT << nl;

    os.writeKeyword("primeUPrimeMean")
        << faItem.primeUPrimeMean_ << token::END_STATEMENT << nl;

    os.writeKeyword("base")
        << faItem.baseTypeNames_[faItem.base_] << token::END_STATEMENT << nl;

    if (faItem.window_ > 0)
    {
        os.writeKeyword("window")
            << faItem.window_ << token::END_STATEMENT << nl;

        if (faItem.windowName_ != "")
        {
            os.writeKeyword("windowName")
                << faItem.windowName_ << token::END_STATEMENT << nl;
        }
    }

    os  << token::END_BLOCK << nl;

    os.check
    (
        "Foam::Ostream& Foam::operator<<"
        "(Foam::Ostream&, const Foam::functionObjects::temporalAveragingItem&)"
    );

    return os;
}


// ************************************************************************* //
