/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

#include "DTFunityLewisEddyDiffusivity.H"
#include "fvmLaplacian.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulenceThermophysicalTransportModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class TurbulenceThermophysicalTransportModel>
DTFunityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>::
DTFunityLewisEddyDiffusivity
(
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo
)
:
    DTFunityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>
    (
        typeName,
        momentumTransport,
        thermo,
        false
    )
{}


template<class TurbulenceThermophysicalTransportModel>
DTFunityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>::
DTFunityLewisEddyDiffusivity
(
    const word& type,
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo,
    const bool allowDefaultPrt
)
:
    unityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>
    (
        type,
        momentumTransport,
        thermo,
        allowDefaultPrt
    ),
    mesh_(momentumTransport.mesh())
{}




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbulenceThermophysicalTransportModels
} // End namespace Foam

// ************************************************************************* //
