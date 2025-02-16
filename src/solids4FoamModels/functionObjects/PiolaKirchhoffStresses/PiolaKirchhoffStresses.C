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

#include "PiolaKirchhoffStresses.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(PiolaKirchhoffStresses, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        PiolaKirchhoffStresses,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::PiolaKirchhoffStresses::writeData()
{
    if 
    (
        runTime_.outputTime() && 
        mesh_.foundObject<volTensorField>("F") &&
        mesh_.foundObject<volSymmTensorField>("sigma")
    )
    {
        // Lookup the deformation gradient and Cachy stress fields
        const volTensorField& F = mesh_.lookupObject<volTensorField>("F");
        const volSymmTensorField& sigma = 
            mesh_.lookupObject<volSymmTensorField>("sigma");

        // Calculate the Jacobian of the deformation gradient
        volScalarField J("J", det(F));

        // Calculate the inverse of the deformation gradient
        volTensorField invF("invF", inv(F));

        // Compute the first (P) and second (S) Piola-Kirchhoff stresses
        volTensorField P("P", J*(inv(F) & sigma));
        volSymmTensorField S("S", symm(P & invF.T()));
        
        // Write the field
        Info<< name_ << ": writing P and S" << endl;
        P.write();
        S.write();
    }
    else
    {
        Info<< name_ << ": F or sigma fields not found!" << endl;
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PiolaKirchhoffStresses::PiolaKirchhoffStresses
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    runTime_(t),
    mesh_
    (
        runTime_.lookupObject<fvMesh>
        (
            dict.lookupOrDefault<word>("region", "region0")
        )
    )
{
    //Info<< "Creating " << this->name() << " function object" << nl
    //    << "    compressionPositive: " << compressionPositive_ << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::PiolaKirchhoffStresses::start()
{
    if (runTime_.outputTime())
    {
        return writeData();
    }

    return true;
}


#if FOAMEXTEND
bool Foam::PiolaKirchhoffStresses::execute(const bool forceWrite)
#else
bool Foam::PiolaKirchhoffStresses::execute()
#endif
{
    if (runTime_.outputTime())
    {
        return writeData();
    }

    return true;
}


bool Foam::PiolaKirchhoffStresses::read(const dictionary& dict)
{
    return true;
}

#ifdef OPENFOAM_NOT_EXTEND
bool Foam::PiolaKirchhoffStresses::write()
{
    return false;
}
#endif

// ************************************************************************* //
