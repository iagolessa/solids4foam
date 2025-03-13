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

Description
    Function that creates and writes principal stretch fields given the 
    deformation gradient tensor.

    Six fields are created:
      - stretchMin, stretchMid and stretchMax: principal stresses ordered from most
        negative to most positive
      - stretchMinDir, stretchMidDir and stretchMaxDir: corresponding principal
        directions
  
Author
    Iago Oliveira, UNESP. 
    Based on code by Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef principalStretchFields_H
#define principalStretchFields_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Set principal stress values and directions
    void calculateEigenValues
    (
        const symmTensor& strainTensor,
        scalar& stretchMax,
        scalar& stretchMid,
        scalar& stretchMin,
        vector& stretchMaxDir,
        vector& stretchMidDir,
        vector& stretchMinDir
    );

    // Calculate and write principal stress value and direction fields
    void writePrincipalStretchFields(const volSymmTensorField& C);


} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
