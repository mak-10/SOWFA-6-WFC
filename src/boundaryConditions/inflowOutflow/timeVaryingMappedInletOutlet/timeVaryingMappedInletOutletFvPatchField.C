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

#include "timeVaryingMappedInletOutletFvPatchField.H"
#include "Time.H"
#include "AverageField.H"
#include "IFstream.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::fileName
Foam::timeVaryingMappedInletOutletFvPatchField<Type>::findFieldFile
(
    const word& timeName
) const
{
    const fileName fieldFileName
    (
        dataDir_/timeName/sampleName_/fieldTableName_
    );

    const fileName typeFieldFileName
    (
        dataDir_/timeName/sampleName_
       /pTraits<Type>::typeName + Field<Type>::typeName
       /fieldTableName_
    );

    if (exists(fieldFileName))
    {
        return fieldFileName;
    }
    else if (exists(typeFieldFileName))
    {
        return typeFieldFileName;
    }
    else
    {
        FatalErrorInFunction
            << "Cannot find field file "
            << fieldFileName << " " << typeFieldFileName
            << exit(FatalError);

        return fileName::null;
    }
}


template<class Type>
void Foam::timeVaryingMappedInletOutletFvPatchField<Type>::checkTable()
{
    // Initialise
    Foam::IOstream::streamFormat fmt = Foam::IOstream::formatEnum(format_);
    if (mapperPtr_.empty())
    {
        // Reread values and interpolate
        const fileName samplePointsFile(dataDir_/pointsName_);

        pointField samplePoints((IFstream(samplePointsFile, fmt)()));

        if (debug)
        {
            Info<< "timeVaryingMappedInletOutletFvPatchField :"
                << " Read " << samplePoints.size() << " sample points from "
                << samplePointsFile << endl;
        }


        // tbd: run-time selection
        bool nearestOnly
        (
           !mapMethod_.empty()
         && mapMethod_ != "planarInterpolation"
        );

        // Allocate the interpolator
        mapperPtr_.reset
        (
            new pointToPointPlanarInterpolation
            (
                samplePoints,
                 this->patch().patch().faceCentres(),
                perturb_,
                nearestOnly
            )
        );

        // Read the times for which data is available
        sampleTimes_ = Time::findTimes(dataDir_);

        if (debug)
        {
            Info<< "timeVaryingMappedInletOutletFvPatchField : In directory "
                << dataDir_ << " found times "
                << pointToPointPlanarInterpolation::timeNames(sampleTimes_)
                << endl;
        }
    }


    // Find current time in sampleTimes
    label lo = -1;
    label hi = -1;

    bool foundTime = mapperPtr_().findTime
    (
        sampleTimes_,
        startSampleTime_,
        this->db().time().value(),
        lo,
        hi
    );

    if (!foundTime)
    {
        FatalErrorInFunction
            << "Cannot find starting sampling values for current time "
            << this->db().time().value() << nl
            << "Have sampling values for times "
            << pointToPointPlanarInterpolation::timeNames(sampleTimes_) << nl
            << "In directory " <<  dataDir_ << " of field " << fieldTableName_
            << exit(FatalError);
    }


    // Update sampled data fields.

    if (lo != startSampleTime_)
    {
        startSampleTime_ = lo;

        if (startSampleTime_ == endSampleTime_)
        {
            // No need to reread since are end values
            if (debug)
            {
                Pout<< "checkTable : Setting startValues to (already read) "
                    << dataDir_/sampleTimes_[startSampleTime_].name()
                    << endl;
            }
            startSampledValues_ = endSampledValues_;
            startAverage_ = endAverage_;
        }
        else
        {
            if (debug)
            {
                Pout<< "checkTable : Reading startValues from "
                    << dataDir_/sampleTimes_[lo].name()
                    << endl;
            }

            // Reread values and interpolate
            const fileName valsFile
            (
                findFieldFile(sampleTimes_[startSampleTime_].name())
            );

            Field<Type> vals;

            if (setAverage_)
            {
                AverageField<Type> avals((IFstream(valsFile, fmt)()));
                vals = avals;
                startAverage_ = avals.average();
            }
            else
            {
                IFstream(valsFile, fmt)() >> vals;
            }

            if (vals.size() != mapperPtr_().sourceSize())
            {
                FatalErrorInFunction
                    << "Number of values (" << vals.size()
                    << ") differs from the number of points ("
                    <<  mapperPtr_().sourceSize()
                    << ") in file " << valsFile << exit(FatalError);
            }

            startSampledValues_ = mapperPtr_().interpolate(vals);
        }
    }

    if (hi != endSampleTime_)
    {
        endSampleTime_ = hi;

        if (endSampleTime_ == -1)
        {
            // endTime no longer valid. Might as well clear endValues.
            if (debug)
            {
                Pout<< "checkTable : Clearing endValues" << endl;
            }
            endSampledValues_.clear();
        }
        else
        {
            if (debug)
            {
                Pout<< "checkTable : Reading endValues from "
                    << dataDir_/sampleTimes_[endSampleTime_].name()
                    << endl;
            }

            // Reread values and interpolate
            const fileName valsFile
            (
                findFieldFile(sampleTimes_[endSampleTime_].name())
            );

            Field<Type> vals;

            if (setAverage_)
            {
                AverageField<Type> avals((IFstream(valsFile, fmt)()));
                vals = avals;
                endAverage_ = avals.average();
            }
            else
            {
                IFstream(valsFile, fmt)() >> vals;
            }

            if (vals.size() != mapperPtr_().sourceSize())
            {
                FatalErrorInFunction
                    << "Number of values (" << vals.size()
                    << ") differs from the number of points ("
                    <<  mapperPtr_().sourceSize()
                    << ") in file " << valsFile << exit(FatalError);
            }

            endSampledValues_ = mapperPtr_().interpolate(vals);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::timeVaryingMappedInletOutletFvPatchField<Type>::
timeVaryingMappedInletOutletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF),
    // inletOutlet:
    phiName_("phi"),
    // timeVaryingMapped:
    fieldTableName_(iF.name()),
    dataDir_(this->db().time().constant()/"boundaryData"/this->patch().name()),
    pointsName_("points"),
    sampleName_(word::null),
    format_("ascii"),
    setAverage_(false),
    fixesValue_(false),
    perturb_(0),
    mapperPtr_(nullptr),
    sampleTimes_(0),
    startSampleTime_(-1),
    startSampledValues_(0),
    startAverage_(Zero),
    endSampleTime_(-1),
    endSampledValues_(0),
    endAverage_(Zero),
    offset_()
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


template<class Type>
Foam::timeVaryingMappedInletOutletFvPatchField<Type>::
timeVaryingMappedInletOutletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    //mixedFvPatchField<Type>(p, iF, dict),
    // skip reading in refGradient
    mixedFvPatchField<Type>(p, iF),
    // inletOutlet:
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    // timeVaryingMapped:
    fieldTableName_(dict.lookupOrDefault("fieldTable", iF.name())),
    dataDir_
    (
        dict.lookupOrDefault
        (
            "dataDir",
            this->db().time().constant()/"boundaryData"/this->patch().name()
        )
    ),
    pointsName_(dict.lookupOrDefault<fileName>("points", "points")),
    sampleName_(dict.lookupOrDefault("sample", word::null)),
    format_(dict.lookupOrDefault<word>("format", "ascii")),
    fixesValue_(dict.lookupOrDefault("fixesValue", false)),
    setAverage_(dict.lookupOrDefault("setAverage", false)),
    perturb_(dict.lookupOrDefault("perturb", 1e-5)),
    mapMethod_
    (
        dict.lookupOrDefault<word>
        (
            "mapMethod",
            "planarInterpolation"
        )
    ),
    mapperPtr_(nullptr),
    sampleTimes_(0),
    startSampleTime_(-1),
    startSampledValues_(0),
    startAverage_(Zero),
    endSampleTime_(-1),
    endSampledValues_(0),
    endAverage_(Zero),
    offset_()
{
    // instead of reading these in from mixedFvPatchField<Type>(p, iF, dict)
    // do the read here so we skip reading refGradient
    this->refValue() = Field<Type>("refValue", dict, p.size());
    this->valueFraction() = scalarField("valueFraction", dict, p.size());

    //-----------------------------------------------------
    //
    // timeVaryingMapped
    //
    dataDir_.expand();
    pointsName_.expand();
    sampleName_.expand();

    if (dict.found("offset"))
    {
        offset_ = Function1<Type>::New("offset", dict);
    }

    if
    (
        mapMethod_ != "planarInterpolation"
     && mapMethod_ != "nearest"
    )
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "mapMethod should be one of 'planarInterpolation'"
            << ", 'nearest'" << exit(FatalIOError);
    }

    if (dict.found("value"))
    {
        fvPatchField<Type>::operator==(Field<Type>("value", dict, p.size()));
    }
    else
    {
        // timeVaryingMapped:
        // Note: we use evaluate() here to trigger updateCoeffs followed
        //       by re-setting of fvatchfield::updated_ flag. This is
        //       so if first use is in the next time step it retriggers
        //       a new update.
        //this->evaluate(Pstream::commsTypes::blocking);

        updateFixedValue(); // update refValues
        fvPatchField<Type>::operator==(this->refValue());
    }

    //-----------------------------------------------------
    //
    // inletOutlet
    //
    // Note: updateCoeffs sets refValue and valueFraction
    //this->refValue() = Zero;
    this->refGrad() = Zero;
    //this->valueFraction() = 0.0;
}


template<class Type>
Foam::timeVaryingMappedInletOutletFvPatchField<Type>::
timeVaryingMappedInletOutletFvPatchField
(
    const timeVaryingMappedInletOutletFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(ptf, p, iF, mapper),
    // inletOutlet:
    phiName_(ptf.phiName_),
    // timeVaryingMapped:
    fieldTableName_(ptf.fieldTableName_),
    dataDir_(ptf.dataDir_),
    pointsName_(ptf.pointsName_),
    sampleName_(ptf.sampleName_),
    format_(ptf.format_),
    setAverage_(ptf.setAverage_),
    fixesValue_(ptf.fixesValue_),
    perturb_(ptf.perturb_),
    mapMethod_(ptf.mapMethod_),
    mapperPtr_(nullptr),
    sampleTimes_(0),
    startSampleTime_(-1),
    startSampledValues_(0),
    startAverage_(Zero),
    endSampleTime_(-1),
    endSampledValues_(0),
    endAverage_(Zero),
    offset_(ptf.offset_, false)
{}


template<class Type>
Foam::timeVaryingMappedInletOutletFvPatchField<Type>::
timeVaryingMappedInletOutletFvPatchField
(
    const timeVaryingMappedInletOutletFvPatchField<Type>& ptf
)
:
    mixedFvPatchField<Type>(ptf),
    // inletOutlet:
    phiName_(ptf.phiName_),
    // timeVaryingMapped:
    fieldTableName_(ptf.fieldTableName_),
    dataDir_(ptf.dataDir_),
    pointsName_(ptf.pointsName_),
    sampleName_(ptf.sampleName_),
    format_(ptf.format_),
    setAverage_(ptf.setAverage_),
    fixesValue_(ptf.fixesValue_),
    perturb_(ptf.perturb_),
    mapMethod_(ptf.mapMethod_),
    mapperPtr_(nullptr),
    sampleTimes_(ptf.sampleTimes_),
    startSampleTime_(ptf.startSampleTime_),
    startSampledValues_(ptf.startSampledValues_),
    startAverage_(ptf.startAverage_),
    endSampleTime_(ptf.endSampleTime_),
    endSampledValues_(ptf.endSampledValues_),
    endAverage_(ptf.endAverage_),
    offset_(ptf.offset_, false)
{}


template<class Type>
Foam::timeVaryingMappedInletOutletFvPatchField<Type>::
timeVaryingMappedInletOutletFvPatchField
(
    const timeVaryingMappedInletOutletFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptf, iF),
    // inletOutlet:
    phiName_(ptf.phiName_),
    // timeVaryingMapped:
    fieldTableName_(ptf.fieldTableName_),
    dataDir_(ptf.dataDir_),
    pointsName_(ptf.pointsName_),
    sampleName_(ptf.sampleName_),
    format_(ptf.format_),
    setAverage_(ptf.setAverage_),
    fixesValue_(ptf.fixesValue_),
    perturb_(ptf.perturb_),
    mapMethod_(ptf.mapMethod_),
    mapperPtr_(nullptr),
    sampleTimes_(ptf.sampleTimes_),
    startSampleTime_(ptf.startSampleTime_),
    startSampledValues_(ptf.startSampledValues_),
    startAverage_(ptf.startAverage_),
    endSampleTime_(ptf.endSampleTime_),
    endSampledValues_(ptf.endSampledValues_),
    endAverage_(ptf.endAverage_),
    offset_(ptf.offset_, false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::timeVaryingMappedInletOutletFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchField<Type>::autoMap(m);
    if (startSampledValues_.size())
    {
        startSampledValues_.autoMap(m);
        endSampledValues_.autoMap(m);
    }
    // Clear interpolator
    mapperPtr_.clear();
    startSampleTime_ = -1;
    endSampleTime_ = -1;
}


template<class Type>
void Foam::timeVaryingMappedInletOutletFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    mixedFvPatchField<Type>::rmap(ptf, addr);

    const timeVaryingMappedInletOutletFvPatchField<Type>& tiptf =
        refCast<const timeVaryingMappedInletOutletFvPatchField<Type>>(ptf);

    startSampledValues_.rmap(tiptf.startSampledValues_, addr);
    endSampledValues_.rmap(tiptf.endSampledValues_, addr);

    // Clear interpolator
    mapperPtr_.clear();
    startSampleTime_ = -1;
    endSampleTime_ = -1;
}


template<class Type>
void Foam::timeVaryingMappedInletOutletFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // The inletOutlet part of the update, which determines the direction of
    // the massflow from phi and updates the valueFraction variable
    updateInletOutlet();

    // The fixedValue part of the update, closely following the updateCoeffs
    // function from timeVaryingMappedFixedValue to update the refValue variable
    updateFixedValue();

    if (debug)
    {
        Info<< "updateCoeffs : " << this->patch().name() << " refValue"
            << " min:" << gMin(this->refValue())
            << " max:" << gMax(this->refValue())
            << " avg:" << gAverage(this->refValue()) << endl;
        Info<< "updateCoeffs : " << this->patch().name() << " valueFraction"
            << " min:" << gMin(this->valueFraction())
            << " max:" << gMax(this->valueFraction())
            << " avg:" << gAverage(this->valueFraction()) << endl;
    }

    // Instead of applying operator==, as is done in the timeVaryingMappedFixedValue,
    // let the field assignment take place in evaluate()
    //
    // Notes:
    // - operator== and operator= are equivalent for fvPatchFields
    // - valueFraction is always 1 or 0--switching between inflow and outflow,
    //   respectively--and is calculated from phi
    // - refGrad == 0 (hard-coded for now), i.e., for outflow, the value is
    //   always extrapolated directly from the interior; this is only used if 
    //   for some reason mixedFvPatchField<Type>::evaluate() is called instead
    //   of timeVaryingMappedInletOutletFvPatchField<Type>::evaluate()
    // - calling updateCoeffs from the parent class sets updated to true
    mixedFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::timeVaryingMappedInletOutletFvPatchField<Type>::updateInletOutlet()
{
    const Field<scalar>& phip =
        this->patch().template lookupPatchField<surfaceScalarField, scalar>
        (
            phiName_
        );

    // Note: this was changed from pos to pos0 between 2.4.x and 6
    //       pos0 returns 1 if positive or zero, otherwise 0
    this->valueFraction() = 1.0 - pos0(phip);
}


template<class Type>
void Foam::timeVaryingMappedInletOutletFvPatchField<Type>::updateFixedValue()
{
    checkTable();

    // Interpolate between the sampled data

    Type wantedAverage;

    if (endSampleTime_ == -1)
    {
        // Only start value
        if (debug)
        {
            Pout<< "updateFixedValue : Sampled, non-interpolated values"
                << " from start time:"
                << sampleTimes_[startSampleTime_].name() << nl;
        }

        //this->operator==(startSampledValues_);
        this->refValue() = startSampledValues_;
        wantedAverage = startAverage_;
    }
    else
    {
        scalar start = sampleTimes_[startSampleTime_].value();
        scalar end = sampleTimes_[endSampleTime_].value();

        scalar s = (this->db().time().value() - start)/(end - start);

        if (debug)
        {
            Pout<< "updateFixedValue : Sampled, interpolated values"
                << " between start time:"
                << sampleTimes_[startSampleTime_].name()
                << " and end time:" << sampleTimes_[endSampleTime_].name()
                << " with weight:" << s << endl;
        }

        //this->operator==((1 - s)*startSampledValues_ + s*endSampledValues_);
        this->refValue() = ((1 - s)*startSampledValues_ + s*endSampledValues_);
        wantedAverage = (1 - s)*startAverage_ + s*endAverage_;
    }

    // Enforce average. Either by scaling (if scaling factor > 0.5) or by
    // offsetting.
    if (setAverage_)
    {
        const Field<Type>& fld = *this;

        Type averagePsi =
            gSum(this->patch().magSf()*fld)
           /gSum(this->patch().magSf());

        if (debug)
        {
            Pout<< "updateFixedValue :"
                << " actual average:" << averagePsi
                << " wanted average:" << wantedAverage
                << endl;
        }

        if (mag(averagePsi) < vSmall)
        {
            // Field too small to scale. Offset instead.
            const Type offset = wantedAverage - averagePsi;
            if (debug)
            {
                Pout<< "updateFixedValue :"
                    << " offsetting with:" << offset << endl;
            }
            this->refValue() += offset;
        }
        else
        {
            const scalar scale = mag(wantedAverage)/mag(averagePsi);

            if (debug)
            {
                Pout<< "updateFixedValue :"
                    << " scaling with:" << scale << endl;
            }
            this->refValue() *= scale;
        }
    }

    // Apply offset to mapped values
    //Info<< "updateFixedValue BEFORE offset: refValue min:" << gMin(this->refValue())
    //    << " max:" << gMax(this->refValue())
    //    << " avg:" << gAverage(this->refValue()) << endl;
    if (offset_.valid())
    {
        const scalar t = this->db().time().timeOutputValue();
        this->refValue() += offset_->value(t);
        //Info<< "updateFixedValue AFTER adding " << offset_->value(t)
        //    << " : refValue min:" << gMin(this->refValue())
        //    << " max:" << gMax(this->refValue())
        //    << " avg:" << gAverage(this->refValue()) << endl;
    }

    if (debug)
    {
        Pout<< "updateFixedValue : set refValue to min:" << gMin(this->refValue())
            << " max:" << gMax(this->refValue())
            << " avg:" << gAverage(this->refValue()) << endl;
    }
}


template<class Type>
void Foam::timeVaryingMappedInletOutletFvPatchField<Type>::evaluate(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    scalar NValueFraction = this->valueFraction().size();
    scalar sumValueFraction = 0;
    forAll(this->valueFraction(), faceI)
    {
        sumValueFraction += this->valueFraction()[faceI];
    }
    reduce(sumValueFraction, sumOp<scalar>());
    reduce(NValueFraction, sumOp<scalar>());
    Info<< "[TVMIO] evaluate() " << this->patch().name()
        << ", avg valueFraction =" << sumValueFraction/NValueFraction
        << endl;

    // This sets the patch to:
    //   valueFraction*refValue
    // + (1 - valueFraction)*(patchInternalField + refGrad*delta)
    Foam::mixedFvPatchField<Type>::evaluate(commsType);
//    fvPatchField<Type>::operator==
//    (
//        this->valueFraction()*this->refValue()
//      +
//        (1.0 - this->valueFraction())*
//        (
//            this->patchInternalField()
//        //+ this->refGrad()/this->patch().deltaCoeffs()
//        )
//    );
//
//    fvPatchField<Type>::evaluate();
}


template<class Type>
void Foam::timeVaryingMappedInletOutletFvPatchField<Type>::write
(
    Ostream& os
) const
{
    fvPatchField<Type>::write(os);
    if (phiName_ != "phi")
    {
        os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
    }

    this->writeEntryIfDifferent
    (
        os,
        "dataDir",
        this->db().time().constant()/"boundaryData"/this->patch().name(),
        dataDir_
    );

    this->writeEntryIfDifferent(os, "points", fileName("points"), pointsName_);

    this->writeEntryIfDifferent(os, "sample", fileName::null, sampleName_);

    this->writeEntryIfDifferent(os, "format", word("ascii"), format_);

    this->writeEntryIfDifferent(os, "setAverage", Switch(false), setAverage_);

    this->writeEntryIfDifferent(os, "fixesValue", Switch(false), fixesValue_);

    this->writeEntryIfDifferent(os, "perturb", scalar(1e-5), perturb_);

    this->writeEntryIfDifferent
    (
        os,
        "fieldTable",
        this->internalField().name(),
        fieldTableName_
    );

    this->writeEntryIfDifferent
    (
        os,
        "mapMethod",
        word("planarInterpolation"),
        mapMethod_
    );

    if (offset_.valid())
    {
        offset_->writeData(os);
    }

    this->refValue().writeEntry("refValue", os);
    //this->refGrad().writeEntry("refGradient", os);
    this->valueFraction().writeEntry("valueFraction", os);

    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::timeVaryingMappedInletOutletFvPatchField<Type>::operator=
(
    const fvPatchField<Type>& ptf
)
{
    // Notes:
    // ------
    // - When bound() is called for a volScalarField, e.g., when solving for
    //   k, boundary values are not updated if this operator= is not defined
    // - bound() is only defined for volScalarField, i.e., a 
    //     GeometricField<scalar, fvPatchField, volMesh>
    //   defined in
    //     finiteVolume/fields/volFields/volFieldsFwd.H
    // - The default GeometricBoundaryField operator= calls the FieldField
    //   operator= and is defined in 
    //     OpenFOAM/fields/GeometricFields/GeometricField/GeometricBoundaryField.C
    //   for which a "FieldField" corresponds to a field of fields, i.e., a
    //   field (list of tensor classes--in this case, boundary patches) of
    //   fields (boundary faces that make up the patch). The operator loops over
    //   fields (boundaries) and performs an assignment, as defined in
    //     OpenFOAM/fields/FieldFields/FieldField/FieldField.C

    //Info<< "[TVMIO] Called operator=" << endl;
    fvPatchField<Type>::operator=(ptf);
}

// ************************************************************************* //
