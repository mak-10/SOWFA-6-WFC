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

#include "windCanopy.H"
#include "volFields.H"
//* #include "/home/mkhan/OpenFOAM/OpenFOAM-6/src/fvOptions/cellSetOption/cellSetOptionI.H" *//

// * * * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::fv::windCanopy::addActuationDiskAxialInertialResistance
(
    vectorField& Usource,
    const labelList& cells,
    const scalarField& Vcells,
    const RhoFieldType& rho,
    const vectorField& U
) const
{
   // scalar a = 1.0 - Cp_/Ct_;
    vector uniDiskDir = diskDir_/mag(diskDir_);
    //tensor E(Zero);
    //E.xx() = uniDiskDir.x();
    //E.yy() = uniDiskDir.y();
    //E.zz() = uniDiskDir.z();

    //vector upU = vector(vGreat, vGreat, vGreat);
    //scalar upRho = vGreat;
    //if (upstreamCellId_ != -1)
    //{
    //    upU =  U[upstreamCellId_];
     //   upRho = rho[upstreamCellId_];
   // }
    //reduce(upU, minOp<vector>());
    //reduce(upRho, minOp<scalar>());


    //scalar T = 2.0*upRho*diskArea_*mag(upU)*a*(1 - a);
    scalar Tt = 0;

    //scalar T = Cft_;

    forAll(cells, i)
    {
	scalar T = 0.5*rho[cells[i]]*sqr(mag(U[cells[i]]))*Cft_*canopyArea_;
	Usource[cells[i]] += ((Vcells[cells[i]]/V())*T*uniDiskDir);
	Tt += (Vcells[cells[i]]/V())*T;
    }
    reduce(Tt, sumOp<scalar>());
    Info << "Total Thrust force in the Canopy " << Tt << nl;
}


// ************************************************************************* //
