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

#include "principalStretches.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "principalStretchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(principalStretches, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        principalStretches,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::principalStretches::writeData()
{
    if (runTime_.outputTime())
    {
        // Lookup stress tensor
        const volTensorField* FPtr = NULL;
        if (mesh_.foundObject<volTensorField>("F"))
        {
            FPtr = &(mesh_.lookupObject<volTensorField>("F"));
        }
        else
        {
            Info<< "No F field" << endl;
        }
        const volTensorField& F = *FPtr;

        // Compute the right Cauch-Green deformation tensor
        volSymmTensorField C("C", symm(F.T() & F));

        // Calculate and write principal stress fields
        writePrincipalStretchFields(C);
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::principalStretches::principalStretches
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
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::principalStretches::start()
{
    if (runTime_.outputTime())
    {
        return writeData();
    }

    return true;
}


#if FOAMEXTEND
bool Foam::principalStretches::execute(const bool forceWrite)
#else
bool Foam::principalStretches::execute()
#endif
{
    if (runTime_.outputTime())
    {
        return writeData();
    }

    return true;
}


bool Foam::principalStretches::read(const dictionary& dict)
{
    return true;
}

#ifdef OPENFOAM_NOT_EXTEND
bool Foam::principalStretches::write()
{
    return false;
}
#endif

// ************************************************************************* //
