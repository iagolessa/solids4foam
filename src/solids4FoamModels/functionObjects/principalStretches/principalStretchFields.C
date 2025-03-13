/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

\*----------------------------------------------------------------------------*/

#include "principalStretchFields.H"

// * * * * * * * * * * * * * * * * * * Functions * * * * * * * * * * * * * * //

void Foam::calculateEigenValues
(
    const symmTensor& strainTensor,
    scalar& stretchMax,
    scalar& stretchMid,
    scalar& stretchMin,
    vector& stretchMaxDir,
    vector& stretchMidDir,
    vector& stretchMinDir
)
{
    const vector eValues = eigenValues(strainTensor);
    const tensor eVectors = eigenVectors(strainTensor);

    label iMax = -1;
    label iMid = -1;
    label iMin = -1;

    // Square root of C eigenvalues are the principal stretches
    const scalar a = Foam::sqrt(eValues[0]);
    const scalar b = Foam::sqrt(eValues[1]);
    const scalar c = Foam::sqrt(eValues[2]);

    if (a < b)
    {
        if (a < c)
        {
            if (b < c)
            {
                // a < b
                // a < c
                // b < c
                // a < b < c
                iMin = 0;
                iMid = 1;
                iMax = 2;
            }
            else
            {
                // a < b
                // a < c
                // b > c
                // a < c < b
                iMin = 0;
                iMid = 2;
                iMax = 1;
            }
        }
        else
        {
            // a < b
            // a > c
            // c < a < b
            iMin = 2;
            iMid = 0;
            iMax = 1;
        }
    }
    else
    {
        if (b < c)
        {
            if (a < c)
            {
                // a > b
                // b < c
                // a < c
                // b < a < c
                iMin = 1;
                iMid = 0;
                iMax = 2;
            }
            else
            {
                // a > b
                // b < c
                // a > c
                // b < c < a
                iMin = 1;
                iMid = 2;
                iMax = 0;
            }
        }
        else
        {
            // a > b
            // b > c
            // c < b < a
            iMin = 2;
            iMid = 1;
            iMax = 0;
        }
    }

    stretchMax = Foam::sqrt(eValues[iMax]);
    stretchMid = Foam::sqrt(eValues[iMid]);
    stretchMin = Foam::sqrt(eValues[iMin]);

    stretchMaxDir = vector(eVectors[3*iMax], eVectors[3*iMax + 1], eVectors[3*iMax + 2]);
    stretchMidDir = vector(eVectors[3*iMid], eVectors[3*iMid + 1], eVectors[3*iMid + 2]);
    stretchMinDir = vector(eVectors[3*iMin], eVectors[3*iMin + 1], eVectors[3*iMin + 2]);
}


void Foam::writePrincipalStretchFields(const volSymmTensorField& C)
{
    const fvMesh& mesh = C.mesh();
    const Time& runTime = mesh.time();

    // Maximum (most positive/tensile) principal stress
    volScalarField stretchMax
    (
        IOobject
        (
            "stretchMax",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("stretchMax", dimless, 0.0)
    );

    volVectorField stretchMaxDir
    (
        IOobject
        (
            "stretchMaxDir",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("stretchMaxDir", dimless, vector::zero)
    );

    // Middle principal stress
    volScalarField stretchMid
    (
        IOobject
        (
            "stretchMid",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("stretchMid", dimless, 0.0)
    );

    volVectorField stretchMidDir
    (
        IOobject
        (
            "stretchMidDir",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("stretchMidDir", dimless, vector::zero)
    );

    // Minimum principal (most negative/compressive) stress
    volScalarField stretchMin
    (
        IOobject
        (
            "stretchMin",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("stretchMin", dimless, 0.0)
    );

    volVectorField stretchMinDir
    (
        IOobject
        (
            "stretchMinDir",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("stretchMinDir", dimless, vector::zero)
    );

    // References to internalFields for efficiency
    const symmTensorField& CI = C.internalField();
#ifdef OPENFOAM_NOT_EXTEND
    scalarField& stretchMaxI = stretchMax.primitiveFieldRef();
    scalarField& stretchMidI = stretchMid.primitiveFieldRef();
    scalarField& stretchMinI = stretchMin.primitiveFieldRef();
    vectorField& stretchMaxDirI = stretchMaxDir.primitiveFieldRef();
    vectorField& stretchMidDirI = stretchMidDir.primitiveFieldRef();
    vectorField& stretchMinDirI = stretchMinDir.primitiveFieldRef();
#else
    scalarField& stretchMaxI = stretchMax.internalField();
    scalarField& stretchMidI = stretchMid.internalField();
    scalarField& stretchMinI = stretchMin.internalField();
    vectorField& stretchMaxDirI = stretchMaxDir.internalField();
    vectorField& stretchMidDirI = stretchMidDir.internalField();
    vectorField& stretchMinDirI = stretchMinDir.internalField();
#endif

    forAll (CI, cellI)
    {
        calculateEigenValues
        (
            CI[cellI],
            stretchMaxI[cellI],
            stretchMidI[cellI],
            stretchMinI[cellI],
            stretchMaxDirI[cellI],
            stretchMidDirI[cellI],
            stretchMinDirI[cellI]
        );
    }

    forAll(C.boundaryField(), patchI)
    {
        if
        (
            !C.boundaryField()[patchI].coupled()
          && mesh.boundaryMesh()[patchI].type() != "empty"
        )
        {
            const symmTensorField& pC = C.boundaryField()[patchI];
#ifdef OPENFOAM_NOT_EXTEND
            scalarField& pStretchMax = stretchMax.boundaryFieldRef()[patchI];
            scalarField& pStretchMid = stretchMid.boundaryFieldRef()[patchI];
            scalarField& pStretchMin = stretchMin.boundaryFieldRef()[patchI];
            vectorField& pStretchMaxDir = stretchMaxDir.boundaryFieldRef()[patchI];
            vectorField& pStretchMidDir = stretchMidDir.boundaryFieldRef()[patchI];
            vectorField& pStretchMinDir = stretchMinDir.boundaryFieldRef()[patchI];
#else
            scalarField& pStretchMax = stretchMax.boundaryField()[patchI];
            scalarField& pStretchMid = stretchMid.boundaryField()[patchI];
            scalarField& pStretchMin = stretchMin.boundaryField()[patchI];
            vectorField& pStretchMaxDir = stretchMaxDir.boundaryField()[patchI];
            vectorField& pStretchMidDir = stretchMidDir.boundaryField()[patchI];
            vectorField& pStretchMinDir = stretchMinDir.boundaryField()[patchI];
#endif

            forAll(pStretchMax, faceI)
            {
                calculateEigenValues
                (
                    pC[faceI],
                    pStretchMax[faceI],
                    pStretchMid[faceI],
                    pStretchMin[faceI],
                    pStretchMaxDir[faceI],
                    pStretchMidDir[faceI],
                    pStretchMinDir[faceI]
                );
            }
        }
    }

    stretchMax.correctBoundaryConditions();
    stretchMid.correctBoundaryConditions();
    stretchMin.correctBoundaryConditions();
    stretchMaxDir.correctBoundaryConditions();
    stretchMidDir.correctBoundaryConditions();
    stretchMinDir.correctBoundaryConditions();

    // Write fields
    Info<< "    Writing stretchMax" << nl
        << "    Writing stretchMid" << nl
        << "    Writing stretchMin" << nl
        << "    Writing stretchMaxDir" << nl
        << "    Writing stretchMidDir" << nl
        << "    Writing stretchMinDir" << endl;

    stretchMax.write();
    stretchMid.write();
    stretchMin.write();
    stretchMaxDir.write();
    stretchMidDir.write();
    stretchMinDir.write();

    Info<< "Principal stretches max = " << gMax(mag(stretchMax)()) << endl;
}


// ************************************************************************* //
