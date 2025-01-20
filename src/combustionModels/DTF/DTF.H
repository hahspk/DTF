/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

Class
    Foam::combustionModels::DTF

Description
    Partially stirred reactor turbulent combustion model.

    This model calculates a finite rate, based on both turbulence and chemistry
    time scales.  Depending on mesh resolution, the Cmix parameter can be used
    to scale the turbulence mixing time scale.

SourceFiles
    DTF.C

\*---------------------------------------------------------------------------*/

#ifndef DTF_H
#define DTF_H

#include "laminar.H"
#include "flameSensor.H"
#include "efficiencyFunction.H"
#include "autoPtr.H"
//#include "basicChemistryModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                            Class DTF Declaration
\*---------------------------------------------------------------------------*/

class DTF
:
    public laminar
{
    // Private Data
        scalar Fmax_;
        volScalarField F_, EF_,tRR_, D_,Fs_;
        autoPtr<flameSensor> flameSensor_; 
        autoPtr<efficiencyFunction> efficiencyFunction_;
        scalar deltaL_;
        scalar SL_;
        scalar upperLimit_;
        scalar lowerLimit_;
        scalar filters_;
    protected:
        // Protected Data
        //- Pointer to chemistry model
        autoPtr<basicChemistryModel> chemistryPtr_;
    public:

    //- Runtime type information
    TypeName("DTF");


    // Constructors

        //- Construct from components
        DTF
        (
            const word& modelType,
            const fluidReactionThermo& thermo,
            const compressibleMomentumTransportModel& turb,
            const word& combustionProperties
        );

        //- Disallow default bitwise copy construction
        DTF(const DTF&);


        //- Destructor
        virtual ~DTF();


    // Member Functions

        //- Correct combustion rate
        virtual void correct();
        
        void flameSpeed();
        //- Fuel consumption rate matrix.
        virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;

        //- Heat release rate [kg/m/s^3]
        virtual tmp<volScalarField> Qdot() const;

        scalar deltaL() const { return deltaL_; }
        scalar SL() const { return SL_; }
        scalar filters() const{return filters_;}
        
        const volScalarField& F() const{return F_;}
        

        //- Update properties from given dictionary
        virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const DTF&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
