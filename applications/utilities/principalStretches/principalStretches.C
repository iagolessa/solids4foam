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

Application
    principalStretches

Description
    Splits up a patch by putting faces in the given bounding box in a new patch

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "volFields.H"
#include "pointFields.H"
#include "principalStretchFields.H"

// * * * * * * * * * * * * *  Main Program * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validOptions.insert("overwrite", "");

#   include "addRegionOption.H"
#   include "addTimeOptions.H"
#   include "setRootCase.H"
#   include "createTime.H"
    // Get times list
    instantList Times = runTime.times();

#   include "checkTimeOptions.H"

    runTime.setTime(Times[startTime], startTime);
    Foam::label endTime = Times.size();

#   include "createNamedMesh.H"

    for (label i = startTime; i < endTime; i++)
    {
        runTime.setTime(Times[i], i);

        Info<< "Time = " << runTime.timeName() << endl;

        mesh.readUpdate();

        // Check if sigma field is found
        IOobject FHeader
        (
            "F",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        if (FHeader.typeHeaderOk<volTensorField>())
        {
            // Read deformation gradient
            const volTensorField F(FHeader, mesh);

            Info << "Computing Cauchy-Green strain tensors\n" << endl;

            // Right Cauchy-Green tensor
            volSymmTensorField C
            (
                IOobject
                (
                    "C",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                 ),
                 symm(F.T() & F)
            );


            Info << "Computing principal stretches\n" << endl;

            // Calculate and write principal stress fields
            writePrincipalStretchFields(C);
        }
        else
        {
            Info<< "    No F field" << endl;
        }

        Info<< endl;
    }

    Info<< "End" << endl;

    return(0);
}


// ************************************************************************* //
