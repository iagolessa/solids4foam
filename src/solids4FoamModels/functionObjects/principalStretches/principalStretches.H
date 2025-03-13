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

Class
    principalStretches

Description
    FunctionObject that outputs the principal stretch fields for nonlinear
    deformation.

    Calculate and write principal stretch fields. It requires the deformation
    gradient tensor F, so it is intended with nonlinear deformation solvers.
    The principal stretches are defined as the eigenvalues of the stretch
    tensor, U, in nonlinear deformations. The stretch tensor is the
    "square-root" of the right Cauchy-Green deformation tensor, C:

        C = U^2

    where C = F.T & F, with F the deformation gradient tensor. Hence, it is
    possible to show that the principal stretches are the square-roots of the
    eingenvalues of C.

    Three vector fields are created: stretchMax, stretchMid, stretchMin.

    @verbatim
    principalStretches1
    {
        // Type of functionObject
        type principalStretches;

        // Where to load it from (if not already in solver)
        // Note: this is not required if the solver already loads this library.
        functionObjectLibs ("libsolids4FoamModels.so");

        // Optional: specifiy mesh region for multi-region solvers, defaults to
        // region0
        region    region0;
    );
    @endverbatim

Author
    Iago Oliveira, UNESP, Brazil. 
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    principalStretches.C

\*---------------------------------------------------------------------------*/

#ifndef principalStretches_H
#define principalStretches_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class principalStretches Declaration
\*---------------------------------------------------------------------------*/

class principalStretches
:
    public functionObject
{
    // Private data

        //- Name
        const word name_;

        //- Reference to main object registry
        const Time& runTime_;

        //- Reference to mesh region
        const fvMesh& mesh_;


    // Private Member Functions

        //- Write data
        bool writeData();

        //- Disallow default bitwise copy construct
        principalStretches
        (
            const principalStretches&
        );

        //- Disallow default bitwise assignment
        void operator=(const principalStretches&);


public:

    //- Runtime type information
    TypeName("principalStretches");


    // Constructors

        //- Construct from components
        principalStretches
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    // Member Functions

        //- start is called at the start of the time-loop
        virtual bool start();

        //- execute is called at each ++ or += of the time-loop
#if FOAMEXTEND
        virtual bool execute(const bool forceWrite);
#else
        virtual bool execute();
#endif

        //- Read and set the function object if its data has changed
        virtual bool read(const dictionary& dict);

#ifdef OPENFOAM_NOT_EXTEND
        //- Write
        virtual bool write();
#endif

#ifndef OPENFOAM_NOT_EXTEND
        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const pointField&)
        {}
#endif
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
