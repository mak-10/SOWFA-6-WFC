/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "drivingForce.H"
#include "interpolate2D.H"
#include "interpolateXY.H"
#include "windRoseToCartesian.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::drivingForce<Type>::updateGivenTimeDepSource_(bool writeIter)
{
    // Interpolate specified source to current time
    Type source = interpolate2D(runTime_.value(),
                                sourceHeightsSpecified_[0],
                                sourceTimesSpecified_,
                                sourceHeightsSpecified_,
                                sourceSpecified_);


    // Apply in all grid cells
    forAll(bodyForce_,cellI)
    {
        bodyForce_[cellI] = source;
    }


    // Write the source information.
    if (writeIter)
    {
        writeSourceHistory_(source);
    }
}


template<class Type>
void Foam::drivingForce<Type>::updateGivenTimeHeightDepSource_(bool writeIter)
{
    // Interpolate specified source values in time and height
    List<Type> source = interpolate2D(runTime_.value(),
                                zPlanes_.planeLocationValues(),
                                sourceTimesSpecified_,
                                sourceHeightsSpecified_,
                                sourceSpecified_);


    // Now go by cell levels and apply the source term
    forAllPlanes(zPlanes_,planeI)
    {
        for (label i = 0; i < zPlanes_.numCellPerPlane()[planeI]; i++)
        {
            label cellI = zPlanes_.planesCellList()[planeI][i];
            bodyForce_[cellI] = source[planeI];
        }
    }


    // Write the column of source information.
    if (writeIter)
    {
        writeSourceHistory_(source);
    }
}


template<class Type>
void Foam::drivingForce<Type>::updateComputedTimeDepSource_(bool writeIter)
{
    // Get the current time step size.
    scalar dt = runTime_.deltaT().value();

    // Interpolate the desired field to current time
    Type desiredField = interpolate2D(runTime_.value(),
                                sourceHeightsSpecified_[0],
                                sourceTimesSpecified_,
                                sourceHeightsSpecified_,
                                sourceSpecified_);
   
    // Convert from speedAndDirection to component if necessary
    Type fldMeanDesired(zeroTensor_());
    if (velocityInputType_ == "component")
    {
        fldMeanDesired = desiredField;
    }
    else if (velocityInputType_ == "speedAndDirection")
    {
        fldMeanDesired = speedDirToComp_(desiredField);
    }


    // Calculate the average field at the closest level to specified
    Type fldMean1 = zPlanes_.average<Type>(field_,hLevels1I);

    // Calculate the average field at the next closest level to specified
    Type fldMean2 = zPlanes_.average<Type>(field_,hLevels2I);

    // Linearly interpolate to get the average field value at the desired height
    Type fldMean = fldMean1 + (((fldMean2 - fldMean1)/(hLevels2 - hLevels1)) * (sourceHeightsSpecified_[0] - hLevels1));
   

    // Compute the source term
    Type ds = (fldMeanDesired - fldMean) / dt;

    // Subtract off any vertical part
    if (setVerticalCompZero_)
    {
        ds = subtractVerticalPart_(ds);
    }

    // Apply the relaxation
    ds *= gain_;

    // Update the source term
    forAll(bodyForce_,cellI)
    {
        bodyForce_[cellI] = ds;
    }

    bodyForce_.correctBoundaryConditions();
    
    // Write the source information
    if (writeIter)
    {
        writeSourceHistory_(ds);
    }
}


template<class Type>
void Foam::drivingForce<Type>::updateComputedTimeHeightDepSource_(bool writeIter)
{
    // Interpolate specified source values in time and height
    List<Type> fldMeanDesired = interpolate2D(runTime_.value(),
                                zPlanes_.planeLocationValues(),
                                sourceTimesSpecified_,
                                sourceHeightsSpecified_,
                                sourceSpecified_);


    // Compute the planar-averaged actual field at each cell level
    List<Type> fldMean = zPlanes_.average<Type>(field_);

    // Compute the error at each cell level
    List<Type> fldError = fldMeanDesired - fldMean;

    // Write the error
    if (writeIter)
    {
        writeErrorHistory_(fldError);
    }

    // Compute the controller action
    List<Type> source = updateController_(fldError);

    // Subtract off any vertical part
    if (setVerticalCompZero_)
    {
        forAllPlanes(zPlanes_,planeI)
        {
            source[planeI] = subtractVerticalPart_(source[planeI]);
        }
    }

    // Blend source terms to constant
    if (blendToConst_)
    {
        blendLinearSlopeDecay_(source);
    }

    // Now go by cell levels and apply the source term
    forAllPlanes(zPlanes_,planeI)
    {
        for (label i = 0; i < zPlanes_.numCellPerPlane()[planeI]; i++)
        {
            label cellI = zPlanes_.planesCellList()[planeI][i];
            bodyForce_[cellI] = source[planeI];
        }
    }
    bodyForce_.correctBoundaryConditions();


    // Write the column of source information.
    if (writeIter)
    {
        writeSourceHistory_(source);
    }
}


template<class Type>
void Foam::drivingForce<Type>::writeSourceHistory_
(
    Type& source
)
{
    // Write the source information.
    if (Pstream::master())
    {
        if (writeSource_)
        {
            if (runTime_.timeIndex() % writeSourceInterval_ == 0)
            {
                sourceHistoryFile_() << runTime_.timeName() << " " << runTime_.deltaT().value() << " " << source << endl;
            }
        }
    }
}


template<class Type>
void Foam::drivingForce<Type>::writeErrorHistory_
(
    Type& error
)
{
    // Write the error information.
    if (Pstream::master())
    {
        if (writeSource_)
        {
            if (runTime_.timeIndex() % writeSourceInterval_ == 0)
            {
                errorHistoryFile_() << runTime_.timeName() << " " << runTime_.deltaT().value() << " " << error << endl;
            }
        }
    }
}


template<class Type>
void Foam::drivingForce<Type>::writeSourceHistory_
(
    List<Type>& source
)
{
    // Write the column of source information.
    if (Pstream::master())
    {
        if (writeSource_)
        {
            if (runTime_.timeIndex() % writeSourceInterval_ == 0)
            {
                sourceHistoryFile_() << runTime_.timeName() << " " << runTime_.deltaT().value();

                forAllPlanes(zPlanes_,planeI)
                {
                   sourceHistoryFile_() << " " << source[planeI];
                }

                sourceHistoryFile_() << endl;
            }
        }
    }
}


template<class Type>
void Foam::drivingForce<Type>::writeErrorHistory_
(
    List<Type>& error
)
{
    // Write the column of source information.
    if (Pstream::master())
    {
        if (writeError_)
        {
            if (runTime_.timeIndex() % writeSourceInterval_ == 0)
            {
                errorHistoryFile_() << runTime_.timeName() << " " << runTime_.deltaT().value();

                forAllPlanes(zPlanes_,planeI)
                {
                   errorHistoryFile_() << " " << error[planeI];
                }

                errorHistoryFile_() << endl;
            }
        }
    }
}


// In general:
// - input type speedAndDirection is not supported
// - subtractVerticalPart does nothing
template<class Type>
Type Foam::drivingForce<Type>::speedDirToComp_
(
    Type desiredField
)
{
    FatalErrorIn
    (
        "Input type 'speedAndDirection' only supported for driving force of type vector"
    ) << abort(FatalError);

    return zeroTensor_();
}


template<class Type>
Type Foam::drivingForce<Type>::subtractVerticalPart_
(
    Type source
)
{
    //Do nothing, just return input
    return source;
}


template<class Type>
void Foam::drivingForce<Type>::blendLinearSlopeDecay_
(
    List<Type>& source
)
{
    Type slope0 = (source[hLevelBlend1] - source[hLevelBlend0])
                / (  zPlanes_.planeLocationValues()[hLevelBlend1]
                   - zPlanes_.planeLocationValues()[hLevelBlend0]);
    Type dslope = -slope0
                / (  zPlanes_.planeLocationValues()[hLevelBlendMax]
                   - zPlanes_.planeLocationValues()[hLevelBlend1]);
    // march from hLevelBlend1 (z >= maxAssimHeight)
    // up to hLevelBlendMax (z >= maxAssimheight + blendThickness)
    // as the slope decreases linearly to 0
    for (label lvl=hLevelBlend1; lvl <= hLevelBlendMax; lvl++)
    {
        Type slope = slope0 + dslope*(  zPlanes_.planeLocationValues()[lvl]
                                      - zPlanes_.planeLocationValues()[hLevelBlend1]);
        source[lvl] = source[lvl-1] + slope*(  zPlanes_.planeLocationValues()[lvl]
                                             - zPlanes_.planeLocationValues()[lvl-1]);
    }
    // set the remaining levels above hLevelBlendMax to the last value
    for (label lvl=hLevelBlendMax+1; lvl < zPlanes_.numberOfPlanes(); lvl++)
    {
        source[lvl] = source[hLevelBlendMax];
    }
}


// Specialization for Type vector
// - input type speedAndDirection is supported
// - subtractVerticalPart is supported
namespace Foam
{
    template<>
    vector Foam::drivingForce<vector>::speedDirToComp_
    (
        vector desiredField
    )
    {
        vector fldMeanDesired = windRoseToCartesian(desiredField.x(),desiredField.y());
        fldMeanDesired.z() = desiredField.z();
        return fldMeanDesired;
    }


    template<>
    vector Foam::drivingForce<vector>::subtractVerticalPart_
    (
        vector source
    )
    {
        vector nUp(vector::zero);
        nUp.z() = 1.0;
        source -= (source & nUp) * nUp;
        return source;
    }
}


template<class Type>
void Foam::drivingForce<Type>::readInputData_(const IOdictionary& ABLProperties)
{
    // PROPERTIES CONCERNING THE SOURCE TERM.
    const dictionary& sourceDict(ABLProperties.subOrEmptyDict(name_ & "Source"));

    // Specify the type of source to use.  The
    // possible types are "given" and "computed".  
    // - The "given" type means that the source values are directly given
    //   and the momentum and temperature fields will react accordingly.  
    // - The "computed" type means that the mean velocity and temperature
    //   are given and the source terms that maintain them are computed. 
    word sourceType(sourceDict.lookup("type"));
    sourceType_ = sourceType;
    

    // If giving the velocity and computing the sources, specify how the velocity
    // is given.  "component" means you enter the x, y, anc z components.
    // "speedAndDirection" means that you enter the horizontal wind speed, horizontal
    // direction, and vertical component.
    if (name_ == "momentum")
    {
        word velocityInputType(sourceDict.lookupOrDefault<word>("inputType","component"));
        velocityInputType_ = velocityInputType;

        bool setVerticalCompZero(sourceDict.lookupOrDefault<bool>("setVerticalCompZero",true));
        setVerticalCompZero_ = setVerticalCompZero;
    }
    else
    {
        velocityInputType_ = "component";
        setVerticalCompZero_ = false;
    }
    

    // Read in the heights at which the sources are given.
    List<scalar> sourceHeightsSpecified = sourceDict.lookup("sourceHeights" & name_);
    sourceHeightsSpecified_ = sourceHeightsSpecified;
    label nSourceHeights = sourceHeightsSpecified_.size();


    // Read in the source table(s) vs. time and height
    readSourceTables_(sourceDict,nSourceHeights);
    

    // Initialize vertical planes
    if (not ((sourceType_ == "given") && (nSourceHeights == 1)) )
    {
        zPlanes_.initialize();
    }

    // Read in the controller gain
    scalar gain(sourceDict.lookupOrDefault<scalar>("gain",1.0));
    gain_ = gain;


    // Profile assimilation
    if ((sourceType_ == "computed") && (nSourceHeights > 1))
    {
        // Write out error profile? 
        bool writeError(sourceDict.lookupOrDefault<bool>("writeError",false));
        writeError_ = writeError;

        // Read in the controller parameters
        scalar alpha(sourceDict.lookupOrDefault<scalar>("alpha",1.0));
        alpha_ = alpha;
    
        scalar timeWindow(sourceDict.lookupOrDefault<scalar>("timeWindow",3600.0));
        timeWindow_ = timeWindow;

        // Transition to constant source terms with height for partial profile assimilation
        blendToConst_ = sourceDict.lookupOrDefault<bool>("blendToConst",false);
    
        if (blendToConst_)
        {
            assimMaxHeight_ = readScalar(sourceDict.lookup("assimMaxHeight"));
            blendThickness_ = sourceDict.lookupOrDefault<scalar>("blendThickness",100.0);
            findBlendLevels_();
        }

        // Smoothing by means of regression curve fitting
        bool regSmoothing(sourceDict.lookupOrDefault<bool>("regSmoothing",true));
        regSmoothing_ = regSmoothing;
    
        if (regSmoothing_)
        {
            label Nreg(sourceDict.lookupOrDefault<label>("regOrder",1));
            Nreg_ = Nreg;
            
            // Read in weights from table
            List<List<scalar> > weightsTable(sourceDict.lookup("weightsTable"));
            // Change from scalar lists to scalar fields
            scalarField heights(weightsTable.size(),0.0);
            scalarField weights(weightsTable.size(),0.0);
            forAll(heights,i)
            {
               heights[i] = weightsTable[i][0];
               weights[i] = weightsTable[i][1];
            }
            // Interpolate to planes
            forAllPlanes(zPlanes_,planeI)
            {
                weights_.append(
                        interpolateXY(zPlanes_.planeLocationValues()[planeI],heights,weights)
                                );
            } 
        }

        // Initialize controller
        initializeController_();
    }


    // If the desired mean wind or temperature is given at only one height, then revert to
    // the old way of specifying the source term.  Find the two grid levels that bracket
    // the given level and interpolate in between them to compute the source term.  So,
    // here, find those two levels.
    if ((sourceType_ == "computed") && (nSourceHeights == 1))
    {
        findSingleForcingHeight_();
    }


    // PROPERTIES CONCERNING GATHERING STATISTICS

    // Gather/write statistics?
    writeSource_ = ABLProperties.lookupOrDefault<bool>("writeSource",true);

    // Statistics gathering/writing frequency?
    writeSourceInterval_ = ABLProperties.lookupOrDefault<int>("writeSourceInterval",1);
}


template<class Type>
void Foam::drivingForce<Type>::readSourceTables_
(
    const dictionary& sourceDict,
    label& nSourceHeights
)
{
    List<List<List<scalar> > > sourceTables;

    for(int i = 0; i < Type::nComponents; i++)
    {
        word sourceTableName = ("sourceTable" & name_) & Type::componentNames[i];

        Info << "Reading " << sourceTableName << endl;
        List<List<scalar> > sourceTable( sourceDict.lookup(sourceTableName) );

        checkSourceTableSize_( sourceTableName, sourceTable, nSourceHeights);

        sourceTables.append(sourceTable);
    }

    // Assuming that the specified times are the same for all components
    label nSourceTimes = sourceTables[0].size();
    Info<< "Processing " << nSourceTimes << " " << name_ << " source times" << endl;
    sourceTimesSpecified_.setSize(nSourceTimes);
    sourceSpecified_ = List<List<Type> >(nSourceTimes, List<Type>(nSourceHeights, Type::zero));

    for(int i = 0; i < nSourceTimes; i++)
    {
        sourceTimesSpecified_[i] = sourceTables[0][i][0];

        for(int j = 0; j < nSourceHeights; j++)
        {
            Type s(Type::zero);

            for(int k = 0; k < Type::nComponents; k++)
            {
                s[k] = sourceTables[k][i][j+1];
            }

            sourceSpecified_[i][j] = s;
        }
    }
}


// Specialization for scalar
// Type scalar does not have nComponent nor componentNames
namespace Foam
{
    template<>
    void drivingForce<scalar>::readSourceTables_
    (
        const dictionary& sourceDict,
        label& nSourceHeights
    )
    {
        word sourceTableName = "sourceTable" & name_;

        Info << "Reading " << sourceTableName << endl;
        List<List<scalar> > sourceTable( sourceDict.lookup(sourceTableName) );

        checkSourceTableSize_( sourceTableName, sourceTable, nSourceHeights);
    
        // Assuming that the specified times are the same for all components
        label nSourceTimes = sourceTable.size();
        Info<< "Processing " << nSourceTimes << " " << name_ << " source times" << endl;
        sourceTimesSpecified_.setSize(nSourceTimes);
        sourceSpecified_ = List<List<scalar> >(nSourceTimes, List<scalar>(nSourceHeights, 0.0));
    
        for(int i = 0; i < nSourceTimes; i++)
        {
            sourceTimesSpecified_[i] = sourceTable[i][0];

            for(int j = 0; j < nSourceHeights; j++)
            {
                sourceSpecified_[i][j] = sourceTable[i][j+1];
            }
        }
    }
}


template<class Type>
void Foam::drivingForce<Type>::checkSourceTableSize_
(
    word& sourceTableName,
    List<List<scalar> >& sourceTable,
    label& nSourceHeights
)
{
    forAll(sourceTable,i)
    {
        if (sourceTable[i].size()-1 != nSourceHeights)
        {
            FatalErrorIn
            (
                "Number of " + sourceTableName + " heights does not match number of heights specified"
            )   << abort(FatalError);
        }
    }
}


template<class Type>
void Foam::drivingForce<Type>::findSingleForcingHeight_()
{
    scalar hLevels1Diff = 0.0;
    scalar hLevels2Diff = 0.0;

    // Find the grid levels closest to the specified height
    List<scalar> hLevelsDiff(zPlanes_.numberOfPlanes());
    forAllPlanes(zPlanes_,planeI)
    {
        hLevelsDiff[planeI] = zPlanes_.planeLocationValues()[planeI] - sourceHeightsSpecified_[0];
    }

    // Find the two levels closest to the specified height
    // Find the closest level
    forAllPlanes(zPlanes_,planeI)
    {
        if (planeI == 0)
        {
            hLevels1I = planeI;
            hLevels1Diff = hLevelsDiff[planeI];
            hLevels1 = zPlanes_.planeLocationValues()[planeI];
        }
        else
        {
            if ( mag(hLevelsDiff[planeI]) < mag(hLevels1Diff) )
            {
                hLevels1I = planeI;
                hLevels1Diff = hLevelsDiff[planeI];
                hLevels1 = zPlanes_.planeLocationValues()[planeI];
            }
        }
    }

    // Find the next closest level
    int j = 0;
    forAllPlanes(zPlanes_,planeI)
    {
        if (planeI != hLevels1I)
        {
            if (j == 0)
            {
                hLevels2I = planeI;
                hLevels2Diff = hLevelsDiff[planeI];
                hLevels2 = zPlanes_.planeLocationValues()[planeI];
            }
            else
            {
                if ( mag(hLevelsDiff[planeI]) < mag(hLevels2Diff) )
                {
                    hLevels2I = planeI;
                    hLevels2Diff = hLevelsDiff[planeI];
                    hLevels2 = zPlanes_.planeLocationValues()[planeI];
                }
            }
            j++;
        }
    }
}


template<class Type>
void Foam::drivingForce<Type>::findBlendLevels_()
{
    hLevelBlend0 = -1;
    hLevelBlend1 = -1;
    hLevelBlendMax = -1;

    // Find the two levels closest to the specified height
    // Find the closest level
    forAllPlanes(zPlanes_,planeI)
    {
        if ( (hLevelBlend1 < 0) &&
             (zPlanes_.planeLocationValues()[planeI] >= assimMaxHeight_) )
        {
            hLevelBlend1 = planeI;
            hLevelBlend0 = planeI - 1;
        }
        else if ( zPlanes_.planeLocationValues()[planeI] >= assimMaxHeight_+blendThickness_ )
        {
            hLevelBlendMax = planeI;
            break;
        }
    }
    Info<< "Blending source terms to constant between z= "
        << zPlanes_.planeLocationValues()[hLevelBlend1] <<  " and "
        << zPlanes_.planeLocationValues()[hLevelBlendMax]
        << endl;
}


template<class Type>
void Foam::drivingForce<Type>::openFiles_()
{
    fileName outputPath(fileName::null);

    // Specify output path
    if (Pstream::parRun())
    {
        outputPath = mesh_.time().path()/".."/"postProcessing"/"sourceHistory"/mesh_.time().timeName();
    }
    else
    {
        outputPath = mesh_.time().path()/"postProcessing"/"sourceHistory"/mesh_.time().timeName();
    }

    if (Pstream::master)
    {
        // Create directory if it does not exist
        mkDir(outputPath);

        word sourceFileName = ("Source" & name_) & "History";
        sourceHistoryFile_.reset(new OFstream(outputPath/sourceFileName));
        
        if (sourceHeightsSpecified_.size() > 1)
        {
            sourceHistoryFile_() << "Heights (m) ";
            forAllPlanes(zPlanes_,planeI)
            {
                sourceHistoryFile_() << zPlanes_.planeLocationValues()[planeI] << " ";
            }
            sourceHistoryFile_() << endl;
        }

        sourceHistoryFile_() << "Time(s)" << " " << "dt (s)" << " " << "source term " << bodyForce_.dimensions() << endl;
        
        //File for saving the error
        if (writeError_)
        {
            word errorFileName = ("Error" & name_) & "History";
            errorHistoryFile_.reset(new OFstream(outputPath/errorFileName));
            
            if (sourceHeightsSpecified_.size() > 1)
            {
                errorHistoryFile_() << "Heights (m) ";
                forAllPlanes(zPlanes_,planeI)
                {
                    errorHistoryFile_() << zPlanes_.planeLocationValues()[planeI] << " ";
                }
                errorHistoryFile_() << endl;
            }

            errorHistoryFile_() << "Time(s)" << " " << "dt (s)" << " " << "error term " << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::drivingForce<Type>::drivingForce
(
    const IOdictionary& dict,
    const word& name,
    const GeometricField<Type, fvPatchField, volMesh>& field
)
:
    // Set name
    name_(name),

    // Set the pointer to runTime
    runTime_(field.time()),

    // Set the pointer to the mesh
    mesh_(field.mesh()),

    // Set the pointer to the velocity field
    field_(field),

    // Compute vertical levels
    zPlanes_(mesh_),

    // Initialize the body force field
    bodyForce_
    (
        IOobject
        (
            "source" & name_,
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensioned<Type>("source" & name_,dimensionSet(field_.dimensions()/dimTime),zeroTensor_())
    ),

    // Initialize height levels and indices for forcing at one height
    hLevels1I(0.0),
    hLevels2I(0.0),
    hLevels1(0.0),
    hLevels2(0.0)
{
    Info << "Creating driving force object for " << name_ << endl;

    readInputData_(dict);
    openFiles_();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::drivingForce<Type>::~drivingForce()
{}


// * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::drivingForce<Type>::update(bool writeIter)
{
    // Source terms are applied directly as given
    if (sourceType_ == "given")
    {
        // Source is only given at one height. In this case,
        // assume the source is set uniformly throughout the domain to the
        // given value as a function of time only.
        if (sourceHeightsSpecified_.size() == 1)
        {
            updateGivenTimeDepSource_(writeIter);
        }

        // Otherwise, set the source as a function of height and time
        else
        {
            updateGivenTimeHeightDepSource_(writeIter);
        }
    }

    // Source terms have to be computed
    else if (sourceType_ == "computed")
    {
        // Field is only specified at one height. In this case,
        // a source term will be computed and applied uniformly in all three
        // dimensions, the same as in the original ABLSolver.
        if (sourceHeightsSpecified_.size() == 1)
        {
            updateComputedTimeDepSource_(writeIter);
        }

        // Otherwise, set the source as a function of height and time
        else
        {
            updateComputedTimeHeightDepSource_(writeIter);
        }
    }
}



// ************************************************************************* //
