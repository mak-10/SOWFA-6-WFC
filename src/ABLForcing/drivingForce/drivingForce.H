/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::drivingForce<Type>

Description
    Determine the driving force

    There are two possible types of source terms:
    - "given": the source values are directly given
      and the flow field will react accordingly.
    - "computed": the mean profiles
      are given and the source terms that maintain them are computed.

    Each object of this class has a name. The name is used to identify
    input parameters. Some special feature are available when name
    equals "momentum"

    General input besides \c type include:
    - \c sourceHeights<name>: heights at which the sources are given
    - \c sourceTable<name>[<componentName>]: source table versus
            time and height (componentwise)

    Additional input when type is \c computed:
    - \c gain: controller gain

    - when name equals "momentum:
        - \c inputType: specify how the velocity is given
            - "component" means the x, y and z components are specified,
            - "speedAndDirection" means horizontal wind speed, cardinal
              direction, and vertical component are specified
        - \c setVerticalCompZero: set vertical part of source term to zero

    - when multiple source heights are specified:
        - \c writeError: Write out error between horizontal averaged profile
          and given input profile
        - \c alpha: Weights of proportional and integral actions are alpha and
          (1-alpha), respectively
        - \c timeWindow: Time window of exponential time filter in integral
          action (s)
        - \c regSmoothing: smooth error by fitting a regression curve
        - when \c regSmoothing is set to true
            - \c regOrder: Order of regression polynomial
            - \c weightsTable: Table with weights versus height
        - \c blendToConst: Blend source terms to a (time-varying) constant
          value above a certain height, to permit assimilation of a partial
          profile; may be applied with direct or indirect profile assimilation
        - when \c blendToConst is set to true
            - \c assimMaxHeight: Height above which the sources are constant
            - \c blendThickness: The vertical distance over which the source
              term transitions to constant with height; the slope changes
              linearly, from the finite difference between the cell at
              assimMaxHeight and just below, to 0 at assimMaxHeight +
              blendThickness

Usage
    \table
        Property                    | Description               | Required  | Default
        type                        | Type of source to use     | yes       |
        sourceHeights<name>         | Source heights            | yes       |
        sourceTable<name>[<...>]    | Source table              | yes       |
        gain                        | Controller gain           | no        | 1.0
        inputType                   | How velocity is specified | no        | "component"
        setVerticalCompZero         | Set vertical part to zero | no        | true
        writeError                  | Write out error profile   | no        | false
        alpha                       | Ratio of P to I action    | no        | 1.0
        timeWindow                  | Time filter window        | no        | 3600.0
        blendToConst                | Constant sources above assimMaxHeight | no | false
        assimMaxHeight              | Constant sources above assimMaxHeight | no | 
        blendThickness              | Vertical distance over which to blend to constant | no |
        regSmoothing                | Regression smoothing      | no        | true
        regOrder                    | Regression order          | no        | 1
        weightsTable                | Table with weights        | no        |
        writeSource                 | Write out source profile  | no        | true
        writeSourceFrequency        | Writing frequency         | no        | 1

SourceFiles
    drivingForce.C

\*---------------------------------------------------------------------------*/

#ifndef drivingForce_H
#define drivingForce_H

#include "fvCFD.H"
#include "IOdictionary.H"
#include "meshPlanes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class drivingForce Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class drivingForce
{
    // Private data


        typedef GeometricField<Type, fvPatchField, volMesh> volFieldType;


        //- Constants
        
            //- Name
            const word name_;

            //- Runtime pointer
            const Time& runTime_;
    
            //- Mesh pointer
            const fvMesh& mesh_;
    
            //- Field pointer
            const volFieldType& field_;


        //- Vertical mesh planes
        meshPlanes zPlanes_;

        //- Driving force field
        volFieldType bodyForce_;

        //- Height level indices for forcing at one height
        label hLevels1I;
        label hLevels2I;

        //- Height levels for forcing at one height
        scalar hLevels1;
        scalar hLevels2;

        //- Height level indices for blending
        label hLevelBlend0;
        label hLevelBlend1;
        label hLevelBlendMax;

        //- Output file pointer
        autoPtr<OFstream> sourceHistoryFile_;

        //- Output file pointer
        autoPtr<OFstream> errorHistoryFile_;

        //- Type of source
        word sourceType_;

        //- Velocity input type (only used when name_ == "momentum")
        word velocityInputType_;

        //- Set vertical component to zero (only for Type vector)
        bool setVerticalCompZero_;

        //- Heights at which the sources are given
        List<scalar> sourceHeightsSpecified_;

        //- Times at which the sources are given
        List<scalar> sourceTimesSpecified_;

        //- Source interpolation table
        List<List<Type> > sourceSpecified_;

        //- Gain of PI controller
        scalar gain_;

        //- Relaxation coefficient (1=P control, 0=I control)
        scalar alpha_;

        //- Time scale of PI controller weighting function
        scalar timeWindow_;

        //- Integrated error
        List<Type> errorInt_;

        //- Direct assimilation of partial profile
        bool blendToConst_;
        scalar assimMaxHeight_;
        scalar blendThickness_;

        //- Weighted regression
        bool regSmoothing_;
        label Nreg_;
        List<scalar> weights_;
        scalarRectangularMatrix Areg_;

        //- Write error profile?
        bool writeError_;

        //- Gather/write source profile?
        bool writeSource_;

        //- Source writing frequency?
        int writeSourceInterval_;


    // Private Member Functions


        //- Return zero object of type for initialisation
        Type zeroTensor_()
        {
            return Type::zero;
        }


        //- Return number of components
        label nComponents_()
        {
            return Type::nComponents;
        }


        //- Read input data from ABLProperties dictionary
        void readInputData_(const IOdictionary& ABLProperties);


        //- Read source tables (componentwise)
        void readSourceTables_
        (
            const dictionary& sourceDict,
            label& nSourceHeights
        );
            

        //- Check size of source table
        void checkSourceTableSize_
        (
            word& sourceTableName,
            List<List<scalar> >& sourceTable,
            label& nSourceHeights
        );

        //- Find grid levels closest to specified height
        void findSingleForcingHeight_();

        //- Find grid levels closest to specified height
        void findBlendLevels_();

        //- Open file streams
        void openFiles_();

        //- Update given time-dependent source term
        void updateGivenTimeDepSource_(bool writeIter);

        //- Update given time-height-dependent source term
        void updateGivenTimeHeightDepSource_(bool writeIter);

        //- Update computed time-dependent source term
        void updateComputedTimeDepSource_(bool writeIter);

        //- Update computed time-height-dependent source term
        void updateComputedTimeHeightDepSource_(bool writeIter);

        //- Initialize controller
        void initializeController_();

        //- Update controller
        List<Type> updateController_(List<Type>& error);

        //- Compute regression coefficients
        scalarRectangularMatrix computeRegressionCoeff_(List<Type>& y);

        //- Construct regression curve
        List<Type> constructRegressionCurve_(scalarRectangularMatrix& beta);

        //- Write source history to file (single value or column of values)
        void writeSourceHistory_(Type& source);
        void writeSourceHistory_(List<Type>& source);

        //- Write error history to file (single value or column of values)
        void writeErrorHistory_(Type& error);
        void writeErrorHistory_(List<Type>& error);

        //- Convert speedAndDirection input to components (only for Type vector)
        Type speedDirToComp_
        (
            Type desiredField
        );

        //- Subtract off vertical part (only for Type vector)
        Type subtractVerticalPart_(Type source);

        //- Blend source terms to constant
        void blendLinearSlopeDecay_(List<Type>& source);

public:

    //- Runtime type information
    ClassName("drivingForce");


    // Constructors
    drivingForce
    (
        const IOdictionary& dict,
        const word& name,
        const volFieldType& field
    );


    // Destructor
    virtual ~drivingForce();


    // Public Member functions

        //- Update source term
        void update(bool writeIter);

        //- Return bodyforce
        volFieldType& force()
        {
            return bodyForce_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// specialization for scalar driving force
template<>
scalar drivingForce<scalar>::zeroTensor_()
{
    return 0.0;
}


template<>
label drivingForce<scalar>::nComponents_()
{
    return 1;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "drivingForce.C"
#   include "drivingForceController.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
