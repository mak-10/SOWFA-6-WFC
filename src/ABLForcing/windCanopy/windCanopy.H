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

Class
    Foam::fv::actuationDiskSource

Description
    Actuation disk source

    Constant values for momentum source for actuation disk
    \f[
        T = 2 \rho A U_{o}^2 a (1-a)
    \f]
    and
    \f[
        U_1 = (1 - a)U_{o}
    \f]

    where:
    \vartable
        A   | disk area
        U_o | upstream velocity
        a   | 1 - Cp/Ct
        U_1 | velocity at the disk
    \endvartable

Usage
    Example usage:
    \verbatim
    fields          (U);        // names of fields to apply source
    diskDir         (-1 0 0);   // disk direction
    Cp              0.1;        // power coefficient
    Ct              0.5;        // thrust coefficient
    diskArea        5.0;        // disk area
    upstreamPoint   (0 0 0);    // upstream point
    \endverbatim


SourceFiles
    windCanopy.C
    windCanopyTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef windCanopy_H
#define windCanopy_H

/*#include "/home/mkhan/OpenFOAM/OpenFOAM-6/src/fvOptions/cellSetOption/cellSetOption.H"*/
#include "cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class actuationDiskSource Declaration
\*---------------------------------------------------------------------------*/

class windCanopy
:
    public cellSetOption
{

protected:

    // Protected data

        //- disk area normal
        vector diskDir_;

        //- Power coefficient
        //scalar Cp_;

        //- Thrust coefficient of the windfarm
        scalar Cft_;

        //- Frontal Area of the wind canopy
        scalar canopyArea_;

        //- Upstream point sample
        point upstreamPoint_;

        //- Upstream cell ID
        label upstreamCellId_;


private:

    // Private Member Functions

        //- Check data
        void checkData() const;

        //- Add resistance to the UEqn
        template<class RhoFieldType>
        void addActuationDiskAxialInertialResistance
        (
            vectorField& Usource,
            const labelList& cells,
            const scalarField& V,
            const RhoFieldType& rho,
            const vectorField& U
        ) const;

        //- Disallow default bitwise copy construct
        windCanopy(const windCanopy&);

        //- Disallow default bitwise assignment
        void operator=(const windCanopy&);


public:

    //- Runtime type information
    TypeName("windCanopy");


    // Constructors

        //- Construct from components
        windCanopy
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~windCanopy()
    {}


    // Member Functions

        // Access

            //- Return Cp
            //scalar Cp() const
            //{
             //   return Cp_;
           // }

            //- Return Cft
            scalar Cft() const
            {
                return Cft_;
            }

            //- Normal disk direction
            const vector& diskDir() const
            {
                return diskDir_;
            }

            //- Frontal Area of the wind canopy
            scalar canopyArea() const
            {
                return canopyArea_;
            }


         // Add explicit and implicit contributions

            //- Source term to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Source term to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldi
            );


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "windCanopyTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
