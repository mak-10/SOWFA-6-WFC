/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "vtkStructuredSetWriter.H"
#include "coordSet.H"
#include "fileName.H"
#include "OFstream.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::vtkStructuredSetWriter<Type>::vtkStructuredSetWriter()
:
    writer<Type>()
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::vtkStructuredSetWriter<Type>::~vtkStructuredSetWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::vtkStructuredSetWriter<Type>::getFileName
(
    const coordSet& points,
    const wordList& valueSetNames
) const
{
    return this->getBaseName(points, valueSetNames) + ".vtk";
}


template<class Type>
void Foam::vtkStructuredSetWriter<Type>::write
(
    const coordSet& points,
    const wordList& valueSetNames,
    const List<const Field<Type>*>& valueSets,
    Ostream& os
) const
{
    /*
    // Work around to not put each file in a separate directory.
    List<word> osComps = os.name().components();
    label osLength = os.name().components().size();
    fileName osNewName;
    for (int i = 0; i < osLength-2; i++)
    {
        osNewName = osNewName + "/" + osComps[i];
    }
    osNewName = osNewName + "/";
    forAll(valueSetNames,i)
    {
        osNewName = osNewName + valueSetNames[i];
        if (i < valueSetNames.size()-1)
        {
            osNewName = osNewName + "_";
        }
    }
    osNewName = osNewName + ".t" + osComps[osComps.size()-2] + ".vtk";

    // Proceed as normal, but with the new file location.    
    OFstream osNew(osNewName); 
    */


    // In the structured format, all we need to know is the min/max 
    // x,y,z coordinates, and how many in each direction, so 
    // we figure that out here.

    List<label> indicesUnique_x;
    List<label> indicesUnique_y;
    List<label> indicesUnique_z;

    List<scalar> points_x(points.size(),0.0);
    List<scalar> points_y(points.size(),0.0);
    List<scalar> points_z(points.size(),0.0);

    forAll(points, i)
    {
        // round to 3 decimal places to help ID unique points
        points_x[i] = std::round(1000.0*points[i].x()) / 1000.0;
        points_y[i] = std::round(1000.0*points[i].y()) / 1000.0;
        points_z[i] = std::round(1000.0*points[i].z()) / 1000.0;
    }

    // Perform sequential stable sorts to get Fortran ordering
    // (x changes fastest, then y, then z)
    List<label> tmpOrder;
    sortedOrder(points_x, tmpOrder);
    List<scalar> new0_points_x(points.size(),0.0);
    List<scalar> new0_points_y(points.size(),0.0);
    List<scalar> new0_points_z(points.size(),0.0);
    List<label> new0_sortOrder(points.size());
    forAll(points, i)
    {
        new0_points_x[i] = points_x[tmpOrder[i]];
        new0_points_y[i] = points_y[tmpOrder[i]];
        new0_points_z[i] = points_z[tmpOrder[i]];
        new0_sortOrder[i] = tmpOrder[i];
        //Info<< "(vtkStructured) x-sorted point " << i
        //    << " (" << new0_points_x[i]
        //    <<  " " << new0_points_y[i]
        //    <<  " " << new0_points_z[i]
        //    << ")" << endl;
    }

    sortedOrder(new0_points_y, tmpOrder);
    List<scalar> new1_points_x(points.size(),0.0);
    List<scalar> new1_points_y(points.size(),0.0);
    List<scalar> new1_points_z(points.size(),0.0);
    List<label> new1_sortOrder(points.size());
    forAll(points, i)
    {
        new1_points_x[i] = new0_points_x[tmpOrder[i]];
        new1_points_y[i] = new0_points_y[tmpOrder[i]];
        new1_points_z[i] = new0_points_z[tmpOrder[i]];
        new1_sortOrder[i] = new0_sortOrder[tmpOrder[i]];
        //Info<< "(vtkStructured) y-sorted point " << i
        //    << " (" << new1_points_x[i]
        //    <<  " " << new1_points_y[i]
        //    <<  " " << new1_points_z[i]
        //    << ")" << endl;
    }

    sortedOrder(new1_points_z, tmpOrder);
    List<label> xyzsort(points.size()); // orig order --> final sorted order
    forAll(points, i)
    {
        points_x[i] = new1_points_x[tmpOrder[i]];
        points_y[i] = new1_points_y[tmpOrder[i]];
        points_z[i] = new1_points_z[tmpOrder[i]];
        xyzsort[i] = new1_sortOrder[tmpOrder[i]];
        //Info<< "(vtkStructured) z-sorted point " << i
        //    << " (" << points_x[i]
        //    <<  " " << points_y[i]
        //    <<  " " << points_z[i]
        //    << ")" << endl;
    }

    uniqueOrder(points_x,indicesUnique_x);
    uniqueOrder(points_y,indicesUnique_y);
    uniqueOrder(points_z,indicesUnique_z);

    List<scalar> pointsUnique_x(indicesUnique_x.size(),0.0);
    forAll(indicesUnique_x, i)
    {
        pointsUnique_x[i] = points_x[indicesUnique_x[i]];
    }

    List<scalar> pointsUnique_y(indicesUnique_y.size(),0.0);
    forAll(indicesUnique_y, i)
    {
        pointsUnique_y[i] = points_y[indicesUnique_y[i]];
    }

    List<scalar> pointsUnique_z(indicesUnique_z.size(),0.0);
    forAll(indicesUnique_z, i)
    {
        pointsUnique_z[i] = points_z[indicesUnique_z[i]];
    }

    // for calculating origin, spacing
    //Info<< "(vtkStructured) pointsUnique_x: " << pointsUnique_x << endl;
    //Info<< "(vtkStructured) pointsUnique_y: " << pointsUnique_y << endl;
    //Info<< "(vtkStructured) pointsUnique_z: " << pointsUnique_z << endl;

    vector minCoords(vector::zero);
    minCoords.x() = min(pointsUnique_x);
    minCoords.y() = min(pointsUnique_y);
    minCoords.z() = min(pointsUnique_z);

    //vector maxCoords(vector::zero);
    //maxCoords.x() = max(pointsUnique_x);
    //maxCoords.y() = max(pointsUnique_y);
    //maxCoords.z() = max(pointsUnique_z);

    label nx = pointsUnique_x.size();
    label ny = pointsUnique_y.size();
    label nz = pointsUnique_z.size();

    scalar dx = 0.0;
    scalar dy = 0.0;
    scalar dz = 0.0;

    if (nx > 1)
    {
        dx = pointsUnique_x[1] - pointsUnique_x[0];

        if (dx < 1e-12)
        {
            FatalErrorIn("vtkStructuredSetWriter<Type>::write(..)")
                << "poinstUnique_x not unique :" << pointsUnique_x
                << exit(FatalError);
        }
    }

    if (ny > 1)
    {
        dy = pointsUnique_y[1] - pointsUnique_y[0];

        if (dy < 1e-12)
        {
            FatalErrorIn("vtkStructuredSetWriter<Type>::write(..)")
                << "poinstUnique_y not unique :" << pointsUnique_y
                << exit(FatalError);
        }
    }

    if (nz > 1)
    {
        dz = pointsUnique_z[1] - pointsUnique_z[0];

        if (dz < 1e-12)
        {
            FatalErrorIn("vtkStructuredSetWriter<Type>::write(..)")
                << "poinstUnique_z not unique :" << pointsUnique_z
                << exit(FatalError);
        }
    }

    if (nx*ny*nz != points.size())
    {
        WarningInFunction
            << "Input points are not structured?" << endl
            << "  nx, ny, nz: " << nx << " " << ny << " " << nz << endl
            << "  number of input points: " << points.size() << endl;
    }

    // Header of the file.
    os.setf(ios_base::fmtflags(ios_base::fixed), ios_base::floatfield);
    os     << "# vtk DataFile Version 3.0" << nl
           << points.name() << nl
           << "ASCII" << nl
           << "DATASET STRUCTURED_POINTS" << nl
           << "DIMENSIONS " << nx << " " << ny << " " << nz << nl
           << "ORIGIN " << minCoords.x() << " " << minCoords.y() << " " << minCoords.z() << nl
           << "SPACING " << dx << " " << dy << " " << dz << nl;

    /*
    forAll(points, i)
    {
        const vector& pt = points[i];
        osNew  << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
    }
    */


    os     << "POINT_DATA " << points.size() << nl
           << " FIELD attributes " << valueSetNames.size() << nl;
        
    forAll(valueSetNames, setI)
    {
     /*
        label nComps = pTraits<Type>::nComponents;
        if (nComps == 1)
        {
            os    << "SCALARS " << "Amb" << " float" << nl;
        }
        else if (nComps == 3)
        {
            os    << "VECTORS " << "Amb" << " float" << nl;
        }
        else if ((nComps == 6) || (nComps == 9))
        {
            os    << "TENSORS " << "Amb" << " float" << nl;
        }
     */
        os << valueSetNames[setI] << ' ' << pTraits<Type>::nComponents << ' '
           << points.size() << " float" << nl;

        const Field<Type>& fld = *valueSets[setI];

        forAll(fld, pointI)
        {
            if (pointI != 0)
            {
                os << ' ';
            }
            
            os.precision(3);
            writer<Type>::write(fld[xyzsort[pointI]], os);
            os << nl;
        }
        os << nl;
    }
}


template<class Type>
void Foam::vtkStructuredSetWriter<Type>::write
(
    const bool writeTracks,
    const PtrList<coordSet>& tracks,
    const wordList& valueSetNames,
    const List<List<Field<Type> > >& valueSets,
    Ostream& os
) const
{
    if (valueSets.size() != valueSetNames.size())
    {
        FatalErrorIn("vtkStructuredSetWriter<Type>::write(..)")
            << "Number of variables:" << valueSetNames.size() << endl
            << "Number of valueSets:" << valueSets.size()
            << exit(FatalError);
    }

    label nTracks = tracks.size();
    label nPoints = 0;
    forAll(tracks, i)
    {
        nPoints += tracks[i].size();
    }

    os  << "# vtk DataFile Version 2.0" << nl
        << tracks[0].name() << nl
        << "ASCII" << nl
        << "DATASET POLYDATA" << nl
        << "POINTS " << nPoints << " float" << nl;

    forAll(tracks, trackI)
    {
        const coordSet& points = tracks[trackI];
        forAll(points, i)
        {
            const vector& pt = points[i];
            os  << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
        }
    }

    if (writeTracks)
    {
        os  << "LINES " << nTracks << ' ' << nPoints+nTracks << nl;

        // Write ids of track points to file
        label globalPtI = 0;
        forAll(tracks, trackI)
        {
            const coordSet& points = tracks[trackI];

            os  << points.size();
            forAll(points, i)
            {
                os  << ' ' << globalPtI;
                globalPtI++;
            }
            os << nl;
        }
    }

    os  << "POINT_DATA " << nPoints << nl
        << " FIELD attributes " << valueSetNames.size() << nl;

    forAll(valueSetNames, setI)
    {
        os  << valueSetNames[setI] << ' ' << pTraits<Type>::nComponents << ' '
            << nPoints << " float" << nl;

        const List<Field<Type> >& fieldVals = valueSets[setI];

        forAll(fieldVals, i)
        {
            const Field<Type>& vals = fieldVals[i];

            forAll(vals, j)
            {
                if (j != 0)
                {
                    os  << ' ';
                }
                writer<Type>::write(vals[j], os);
            }
            os  << nl;
        }
    }
}


// ************************************************************************* //
