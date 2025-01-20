/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2022 OpenFOAM Foundation
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

#include "DTFunityLewisFourier.H"
#include "fvmLaplacian.H"
#include "fvcSnGrad.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarThermophysicalTransportModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicThermophysicalTransportModel>
DTFunityLewisFourier<BasicThermophysicalTransportModel>::DTFunityLewisFourier
(
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo
)
:
    DTFunityLewisFourier
    (
        typeName,
        momentumTransport,
        thermo
    )
{}


template<class BasicThermophysicalTransportModel>
DTFunityLewisFourier<BasicThermophysicalTransportModel>::DTFunityLewisFourier
(
    const word& type,
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo
)
:
    unityLewisFourier<BasicThermophysicalTransportModel>
    (
        type,
        momentumTransport,
        thermo
    ),
    mesh_(momentumTransport.mesh())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



template<class BasicThermophysicalTransportModel>
tmp<surfaceScalarField>
DTFunityLewisFourier<BasicThermophysicalTransportModel>::q() const
{
const volScalarField& EF = mesh_.lookupObject<volScalarField>("EF");
    return surfaceScalarField::New
    (
        IOobject::groupName
        (
            "q",
            this->momentumTransport().alphaRhoPhi().group()
        ),
       -fvc::interpolate(this->alpha()*this->thermo().alphahe()*EF)
       *fvc::snGrad(this->thermo().he())
    );
}


template<class BasicThermophysicalTransportModel>
tmp<fvScalarMatrix>
DTFunityLewisFourier<BasicThermophysicalTransportModel>::
divq(volScalarField& he) const
{
const volScalarField& EF = mesh_.lookupObject<volScalarField>("EF");
    return -fvm::laplacian(this->alpha()*this->thermo().alphahe()*EF, he);
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace laminarThermophysicalTransportModels
} // End namespace Foam

// ************************************************************************* //
