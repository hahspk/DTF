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

#include "efficiencyFunction.H"
#include "DTF.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(efficiencyFunction, 0);

    defineRunTimeSelectionTable(efficiencyFunction, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::efficiencyFunction::efficiencyFunction
(
    const dictionary& dict,
    const fvMesh& mesh,
    const compressibleMomentumTransportModel& turb,
    scalar Fmax,
    const Foam::combustionModels::DTF& dtfModel
)
:
    mesh_(mesh),
    E_(
         IOobject
         (
             "E",
             //mesh.time().name(),
             mesh().time().timeName(),
             mesh(),
             IOobject::NO_READ,
             IOobject::NO_WRITE
         ),
         mesh,
         dimensionedScalar(dimless, 1)
     ),
     F_(Fmax),
     turb_(turb),
     dtfModel_(dtfModel)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::efficiencyFunction::~efficiencyFunction()
{}


// ************************************************************************* //
