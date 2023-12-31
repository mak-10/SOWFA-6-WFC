    // Re-read blending information
    List<List<scalar> > profileTable(mesh.schemesDict().subDict("schemeBlending").lookup("faceSizeBlendingTable"));
    forAll(faceAreaList,i)
    {
        faceAreaList[i] = profileTable[i][0];
        UBlendingList[i] = profileTable[i][1];
        TBlendingList[i] = profileTable[i][2];
    }

    zBlending1 = readScalar(mesh.schemesDict().subDict("schemeBlending").lookup("heightBlending_z1"));
    zBlending2 = readScalar(mesh.schemesDict().subDict("schemeBlending").lookup("heightBlending_z2"));

    blendingFactorUz1 = readScalar(mesh.schemesDict().subDict("schemeBlending").lookup("heightBlendingFactorU_z1"));
    blendingFactorTz1 = readScalar(mesh.schemesDict().subDict("schemeBlending").lookup("heightBlendingFactorT_z1"));
    blendingFactorUz2 = readScalar(mesh.schemesDict().subDict("schemeBlending").lookup("heightBlendingFactorU_z2"));
    blendingFactorTz2 = readScalar(mesh.schemesDict().subDict("schemeBlending").lookup("heightBlendingFactorT_z2"));


    bool useWallDistZ(mesh.schemesDict().subDict("schemeBlending").lookupOrDefault<bool>("useWallDistZ",false));



    // Check to see if the blending factors or face sizes have changed.  If so, recalculate the U and T
    // blending factor field.
    if ((profileTable != profileTableOld) ||
        (zBlending1 != zBlending1Old) ||
        (zBlending2 != zBlending2Old) ||
        (blendingFactorUz1 != blendingFactorUz1Old) ||
        (blendingFactorTz1 != blendingFactorTz1Old) ||
        (blendingFactorUz2 != blendingFactorUz2Old) ||
        (blendingFactorTz2 != blendingFactorTz2Old) ||
        (useWallDistZ != useWallDistZOld))

    {


        // Write a message to the log file.
        Info << "Updating scheme blending factor field" << endl;


        // Update the old value of these quantities for checking for updated files.
        profileTableOld = profileTable;

        zBlending1Old = zBlending1;
        zBlending2Old = zBlending2;

        blendingFactorUz1Old = blendingFactorUz1;
        blendingFactorTz1Old = blendingFactorTz1;
        blendingFactorUz2Old = blendingFactorUz2;
        blendingFactorTz2Old = blendingFactorTz2;

        useWallDistZOld = useWallDistZ;

        // Get distance from the wall
        vector up = vector::zero;
        up.z() = 1.0;
        surfaceScalarField dFace = mesh.Cf() & up;
        if (useWallDistZ)
        {
            Info << "Calculating wall distance" << endl;
            wallDist d(mesh);
            dFace = fvc::interpolate(d.y());
        }





        // Update the blending factor fields.
        // internal field
        forAll(UBlendingFactor, faceI)
        {
            scalar area = mesh.magSf()[faceI];

            UBlendingFactor[faceI] = max(min(interpolateXY(area,faceAreaList,UBlendingList),1.0),0.0);
            TBlendingFactor[faceI] = max(min(interpolateXY(area,faceAreaList,TBlendingList),1.0),0.0);




            scalar z = 0.0;
            if (useWallDistZ)
            {
                z = dFace[faceI];
            }
            else
            {
                z = mesh.Cf()[faceI].z();
            }
  
            if ((z > zBlending1) && (z < zBlending2))
            {
                UBlendingFactor[faceI] = blendingFactorUz2 +
                                         0.5 * (blendingFactorUz1 - blendingFactorUz2) *
                                        (1.0 + Foam::cos(((z - zBlending1)/(zBlending2 - zBlending1))*Foam::constant::mathematical::pi));
                TBlendingFactor[faceI] = blendingFactorTz2 +
                                         0.5 * (blendingFactorTz1 - blendingFactorTz2) *
                                        (1.0 + Foam::cos(((z - zBlending1)/(zBlending2 - zBlending1))*Foam::constant::mathematical::pi));
            }
            else if ( z > zBlending2 )
            {
                UBlendingFactor[faceI] = blendingFactorUz2;
                TBlendingFactor[faceI] = blendingFactorTz2;
            }

        }

        // boundary field
        forAll(UBlendingFactor.boundaryFieldRef(), patchI)
        {
            forAll(UBlendingFactor.boundaryFieldRef()[patchI], faceI)
            {
                scalar area = mesh.boundary()[patchI].magSf()[faceI];

                UBlendingFactor.boundaryFieldRef()[patchI][faceI] = max(min(interpolateXY(area,faceAreaList,UBlendingList),1.0),0.0);
                TBlendingFactor.boundaryFieldRef()[patchI][faceI] = max(min(interpolateXY(area,faceAreaList,TBlendingList),1.0),0.0);

                scalar z = 0.0;
                if (useWallDistZ)
                {
                    z = dFace.boundaryField()[patchI][faceI];
                }
                else
                {
                    z = mesh.boundary()[patchI].Cf()[faceI].z();
                }


                if ((z > zBlending1) && (z < zBlending2))
                {
                    UBlendingFactor.boundaryFieldRef()[patchI][faceI] = blendingFactorUz2 +
                                             0.5 * (blendingFactorUz1 - blendingFactorUz2) *
                                            (1.0 + Foam::cos(((z - zBlending1)/(zBlending2 - zBlending1))*Foam::constant::mathematical::pi));
                    TBlendingFactor.boundaryFieldRef()[patchI][faceI] = blendingFactorTz2 +
                                             0.5 * (blendingFactorTz1 - blendingFactorTz2) *
                                            (1.0 + Foam::cos(((z - zBlending1)/(zBlending2 - zBlending1))*Foam::constant::mathematical::pi));
                }
                else if ( z > zBlending2 )
                {
                    UBlendingFactor.boundaryFieldRef()[patchI][faceI] = blendingFactorUz2;
                    TBlendingFactor.boundaryFieldRef()[patchI][faceI] = blendingFactorTz2;
                }

            }
        }
    }
