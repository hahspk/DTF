/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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

#include "tanhSensor.H"
#include "addToRunTimeSelectionTable.H"
#include "combustionModel.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace flameSensorModels
{
    defineTypeNameAndDebug(tanhSensor, 0);

    addToRunTimeSelectionTable
    (
        flameSensor,
        tanhSensor,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flameSensorModels::tanhSensor::tanhSensor
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    flameSensor(dict, mesh),
    coeffsDict_(dict.subDict("tanhCoeffs")),
    beta_(coeffsDict_.lookup<scalar>("beta")),
    n_filters_(coeffsDict_.lookup<scalar>("n_filters")),
    sFilter_(mesh_)
{

}

void Foam::flameSensorModels::tanhSensor::correct() {
const word modelName
(
 IOobject::groupName
 (
     combustionModel::combustionPropertiesName,
     ""
 )
);
tmp<volScalarField> tqdot = mesh_.lookupObject<combustionModel>(modelName).Qdot();
/////////////////////////////////simpleFilter/////////////////////////////////
    for(int i = 0; i < n_filters_; i++) {
      tqdot = sFilter_(tqdot);
     } 
volScalarField qdot = Foam::mag(tqdot.ref());

/*
forAll(qdot, cellI)
{
    qdot[cellI] = max(qdot[cellI], 100);
}
*/

dimensionedScalar maxValue(qdot.dimensions(),gMax(qdot));
Info << "qdot maximum: "<< gMax(qdot) << "Minimum : "<< gMin(qdot)<< endl;

S_ = tanh(beta_*qdot/maxValue);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::flameSensorModels::tanhSensor::~tanhSensor()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// ************************************************************************* //
