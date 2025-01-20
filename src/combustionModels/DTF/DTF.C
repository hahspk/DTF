/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
 
#include "DTF.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{
    defineTypeNameAndDebug(DTF, 0);
    addToRunTimeSelectionTable(combustionModel, DTF, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModels::DTF::DTF
(
    const word& modelType,
    const fluidReactionThermo& thermo,
    const compressibleMomentumTransportModel& turb,
    const word& combustionProperties
)
:
    laminar(modelType, thermo, turb, combustionProperties),
    Fmax_(this->coeffs().template lookup<scalar>("Fmax")),
    F_
     (
         IOobject
         (
             "F",
             //this->mesh().time().name(),
             this->mesh().time().timeName(),
             this->mesh(),
             IOobject::NO_READ,
             IOobject::NO_WRITE
         ),
         this->mesh(),
         dimensionedScalar(dimless, 1)
     ),
    EF_
     (
         IOobject
         (
             "EF",
             //this->mesh().time().name(),
             this->mesh().time().timeName(),
             this->mesh(),
             IOobject::NO_READ,
             IOobject::NO_WRITE
         ),
         this->mesh(),
         dimensionedScalar(dimless, 1)
     ), 
    tRR_
    (
        IOobject
        (
            "tRR",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("tRR", dimMass/dimVolume/dimTime, Zero)
    ),
    D_
    (
        IOobject
        (
            "D",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("D", dimLength*dimLength/dimTime, Zero)
    ),
       Fs_
     (
         IOobject
         (
             "Fs",
             //this->mesh().time().name(),
             this->mesh().time().timeName(),
             this->mesh(),
             IOobject::NO_READ,
             IOobject::NO_WRITE
         ),
         this->mesh(),
         dimensionedScalar(dimless, 1)
     ),
    flameSensor_(flameSensor::New(coeffs(),turb.mesh())),
    efficiencyFunction_(efficiencyFunction::New(turb.mesh(),coeffs(),turb,Fmax_, *this)),
    deltaL_(0.001),
    SL_(3),
    upperLimit_(this->coeffs().template lookup<scalar>("upperLimit")),
    lowerLimit_(this->coeffs().template lookup<scalar>("lowerLimit")),
    filters_(this->coeffs().template lookup<scalar>("filters")),
    chemistryPtr_(basicChemistryModel::New(thermo))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::combustionModels::DTF::~DTF()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::combustionModels::DTF::correct()
{
    laminar::correct();
    flameSensor_->correct();
    efficiencyFunction_->correct();
    ////////////////////////////////////////////////////////////////////////// 
    ///FlameSPeed Switch, by default it is true. Set to False incase of cold (non-reactive Flows)
    const bool FlameSpeed = this->coeffs().lookupOrDefault("FlameSpeed", true); 
    if (FlameSpeed)
    {
        flameSpeed();
    }
    /////////////////////////Cell_Size_Implementation/////////////////////////
    //loop will calculate the dyanmic thickened factor according to the cell size and flame thickness.
    forAll(F_, celli) {
    //Cell Size
    scalar delta_g = pow(mesh_.V()[celli], 1.0 / static_cast<scalar>(mesh_.nGeometricD()));
    if (deltaL_ > SMALL) {
        //thickening Factor
        scalar Fvalue = (delta_g* this->coeffs().template lookup<scalar>("N") / (deltaL_ ));
        //switch to choose the minimum of thickening factor, between calculated thickening factor and user defined maximum value.
        Fs_[celli] = min(Fvalue, Fmax_);
        
        scalar sensorValue = flameSensor_->S()[celli];
        //Applying thickening to the flame
        F_[celli] = 1 + (Fs_[celli] - 1) * sensorValue;
        //check if flame thickening is less then 1.
        if (F_[celli] < 1) {
            F_[celli] = 1;
        }
        
        /////////////////////////Log/////////////////////////////////////
        Info << "Maximum Thickening Factor, F_: " << Foam::gMax(tRR_) << endl;
        Info << "Minimum Thickening Factor, F_: " << Foam::gMin(tRR_) << endl;

        } else {
            F_[celli] = 1.0;
        ////////////////////////////log//////////////////////////////////////////
        if (celli % 100 == 0)  // Progress indicator and debugging every 100 cells
        {
            Info << "Cell " << celli << ": deltaL_ too small, set F_[" << celli << "] to " << "Fmin" << endl;
        }
        }
    }
        Info << "Maximum Thickening Factor, Fs_: " << Foam::gMax(Fs_) << endl;
        Info << "Minimum Thickening Factor, Fs_: " << Foam::gMin(Fs_) << endl;
        Info << "Maximum Thickening Factor, F_: " << Foam::gMax(F_) << endl;
        Info << "Minimum Thickening Factor, F_: " << Foam::gMin(F_) << endl;
    //////////////////////////Efficiency into Thickening Factor///////////
    forAll(EF_, celli) {
    EF_[celli] = efficiencyFunction_->E()[celli] * F_[celli];
    }
        Info << "Maximum Efficiency Factor, E_: " << Foam::gMax(efficiencyFunction_->E()) << endl;
        Info << "Minimum Efficiency Factor, E_: " << Foam::gMin(efficiencyFunction_->E()) << endl;
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::DTF::R(volScalarField& Y) const
{
    return efficiencyFunction_->E()/F_*laminar::R(Y);
}

Foam::tmp<Foam::volScalarField>
Foam::combustionModels::DTF::Qdot() const
{
    return volScalarField::New
    (
        this->thermo().phasePropertyName(typeName + ":Qdot"),
        efficiencyFunction_->E()/F_*laminar::Qdot()
    );
}


void Foam::combustionModels::DTF::flameSpeed()
{
    // get number of species 
    const label nSpecie = chemistryPtr_->nSpecie();
    Info << "Number of species: " << nSpecie << endl;

    //get number of reactions
    const label nReaction = chemistryPtr_->nReaction();
    Info << "Number of reactions: " << nReaction << endl;

    ////////////////////////// Reaction Rate Calculations /////////////////////
    
    // Reset reactionRate for all cells to zerro
    forAll(tRR_.internalField(), celli)
    {
        tRR_[celli] = 0.0;
    }

    // calculate the reaction rate for all the recations and species
    for (label ri = 0; ri < nReaction; ++ri)
    {
        for (label si=0; si<nSpecie; si++)
        {
        volScalarField::Internal RR = Foam::mag(chemistryPtr_->calculateRR(ri, si));
        forAll(RR, celli)
        {
            if (!std::isnan(F_[celli]) && F_[celli] > SMALL)
            {
                tRR_[celli] +=  RR[celli];
            }
        }
        }
    }
    volScalarField smoothed_tRR = tRR_; // Initialize smoothed field with raw tRR_

    //smoothing filter
    fvc::smooth(smoothed_tRR, filters_); 
    tRR_ = smoothed_tRR;  
    ///////////////////////// Log //////////////////////////
    dimensionedScalar maxSumRR_(tRR_.dimensions(),gMax(tRR_)); 
    Info << "Maximum reaction rate, maxSumRR_: " << maxSumRR_ << endl;
    Info << "Minimum reaction rate, minSumRR_: " << Foam::gMin(tRR_) << endl;
    //////////////////////////////////////////////////////////////////////////////
    //Flame Sensor
    volScalarField q_sensor = flameSensor_->S();
    scalar sumReactionRate = 0.0;
    scalar count = 0.0;
    Info << "Internal field size: " << q_sensor.internalField().size() << endl;
    //Extracting values reaction rates between the set values of flame sensor
    forAll(q_sensor.internalField(), celli) 
    {
        if (q_sensor[celli] >= lowerLimit_ && q_sensor[celli] <= upperLimit_)
        {
            // Summing reaction rates and counting cells
            sumReactionRate += tRR_[celli];
            count++;
        }
    }
    // Parallel reduction for sum and count
    scalar globalSumReactionRate = sumReactionRate;
    scalar globalCount = count;
    reduce(globalSumReactionRate, sumOp<scalar>());
    reduce(globalCount, sumOp<scalar>());

    // Compute average reaction rate across all processors
    scalar avgReactionRate = (globalCount > 0) ? (globalSumReactionRate / globalCount) : 0.0;

    Info << "Average reaction rate (based on sensor): " << avgReactionRate << endl;

    ///////////////////////// Diffusion Calculations /////////////////////////
    const volScalarField& rhoField = Foam::combustionModel::thermo_.rho();
    const volScalarField& kappaField = Foam::combustionModel::thermo_.kappa();
    const volScalarField& CpField = Foam::combustionModel::thermo_.Cp();

    forAll(this->mesh().V(), celli)
    {
        scalar rho = rhoField[celli];
        scalar lambda = kappaField[celli];
        scalar cp = CpField[celli];

        if (rho > SMALL && cp > SMALL)
        {
            D_[celli] =  lambda / ((rho ) * (cp));   
        }
        else
        {
            D_[celli] = 0;
        }
    }
        Info << "Maximum Diffusion, D_:" << Foam::gMax(D_) << endl;
        Info << "Minimum Diffusion, D_:" << Foam::gMin(D_) << endl;

        //Smootinh Filter
        volScalarField smoothed_D = D_; 
        fvc::smooth(smoothed_D, filters_); 
        D_ = smoothed_D;  
        scalar maxD = Foam::gMax(D_);
        Info << "Maximum Diffusion: smoothed_D" << maxD << endl;
        Info << "Minimum Diffusion: smoothed_D" << Foam::gMin(D_) << endl;

        scalar sumD = 0.0;
        scalar countD = 0.0; 
        //Extracting Values of diffusion between the set values of flame sensor
        forAll(q_sensor.internalField(), celli)
        {
            if (q_sensor[celli] >= lowerLimit_ && q_sensor[celli] <= upperLimit_)
            {
                sumD += D_[celli];
                countD++;
            }
        }

        // Parallel reduction for sumD and countD
        scalar globalSumD = sumD;
        scalar globalCountD = countD;
        reduce(globalSumD, sumOp<scalar>());
        reduce(globalCountD, sumOp<scalar>());

        // Compute average D across all processors
        scalar avgD = (globalCountD > 0) ? (globalSumD / globalCountD) : 0.0;

        Info << "Average Diffusion (based on sensor): " << avgD << endl;

        /////////////////////// SL_ and deltaL_ Calculation ///////////////////////
        Info << "Starting Lamianr flame speed (SL_) and Laminar flame thickness (deltaL_) calculation..." << endl;

            if ( avgReactionRate > SMALL && avgD > SMALL )
            {
                SL_ = sqrt(avgReactionRate * avgD);
                Info << "Laminar Flame Speed: " << SL_ << endl;
            
                deltaL_ = avgD /SL_;

                Info << "Laminar Flame thickness: " << SL_ << endl;

                if (SL_ > 50)
                {
                    Info <<  " SL_ value " << SL_ << " exceeds 50, clamping to 50." << endl;
                    SL_ = 50;
                }

                if (deltaL_ > 10)
                {
                    Info << "deltaL_ value " << deltaL_ << " exceeds 10, clamping to 10." << endl;
                    deltaL_ = 10;
                }
            }
            else
            {
                SL_ = this->coeffs().template lookup<scalar>("SL");
                deltaL_ = this->coeffs().template lookup<scalar>("deltaL");
            }
        Info << "Completed Lamianr flame speed  (SL_) and Laminar flame thickness (deltaL_) calculation." << endl;
}


bool Foam::combustionModels::DTF::read()
{
    if (laminar::read())
    {
        this->coeffs().lookup("Fmax") >> Fmax_;
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
