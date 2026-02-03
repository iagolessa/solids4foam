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
    principalStresses

Description
    Splits up a patch by putting faces in the given bounding box in a new patch

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "volFields.H"
#include "pointFields.H"
#include "principalStressFields.H"

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

    //- Flag to signify if compression is considered positive
    const Switch compressionPositive = false;

    for (label i = startTime; i < endTime; i++)
    {
        runTime.setTime(Times[i], i);

        Info<< "Time = " << runTime.timeName() << endl;

        mesh.readUpdate();

        // Check if sigma field is found
        IOobject sigmaHeader
        (
            "sigma",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        if (sigmaHeader.typeHeaderOk<volSymmTensorField>())
        {
            const volSymmTensorField sigma(sigmaHeader, mesh);

            Info << "Computing principal stresses\n" << endl;

            // Calculate and write principal stress fields
            writePrincipalStressFields(sigma, compressionPositive);
        }
        else
        {
            Info<< "    No sigma field" << endl;
        }

        Info<< endl;
    }

    Info<< "End" << endl;

    return(0);

    Info<< nl << "End" << nl << endl;

    return(0);
}


// ************************************************************************* //
