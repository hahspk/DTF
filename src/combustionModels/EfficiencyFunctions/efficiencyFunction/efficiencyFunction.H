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

Namespace
    Foam::efficiencyFunctionModels

Description
    Namespace for laminar flame speed models


Class
    Foam::efficiencyFunction

Description
    Abstract class for laminar flame speed

SourceFiles
    efficiencyFunction.C
    efficiencyFunctionNew.C

\*---------------------------------------------------------------------------*/

#ifndef efficiencyFunction_H
#define efficiencyFunction_H

#include "dictionary.H"
#include "volFields.H"
#include "fvMesh.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "compressibleMomentumTransportModels.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam {
namespace combustionModels {
    class DTF;
}
} 

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class efficiencyFunction Declaration
\*---------------------------------------------------------------------------*/



class efficiencyFunction
{
protected:

    // Protected data

        const fvMesh& mesh_;
        volScalarField E_;
        scalar F_;
        const compressibleMomentumTransportModel& turb_;
        const Foam::combustionModels::DTF& dtfModel_;


public:

    //- Runtime type information
    TypeName("efficiencyFunction");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            efficiencyFunction,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& mesh,
                const compressibleMomentumTransportModel& turb,
                scalar F,
                const Foam::combustionModels::DTF& dtfModel
            ),
            (dict, mesh, turb, F, dtfModel)
        );


    // Constructors

        //- Construct from dictionary and psiuReactionThermo
        efficiencyFunction
        (
            const dictionary&,
            const fvMesh& mesh,
            const compressibleMomentumTransportModel& turb,
            scalar F,
            const Foam::combustionModels::DTF& dtfModel
        );

        //- Disallow default bitwise copy construction
        efficiencyFunction(const efficiencyFunction&) = delete;


    // Selector

        static autoPtr<efficiencyFunction> New
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const compressibleMomentumTransportModel& turb,
            scalar F,
            const Foam::combustionModels::DTF& dtfModel
        );

    // Member functions
    
    virtual void correct() = 0;
    
    const volScalarField& E() const{
    	return E_;
    }
    
    //- Destructor
    virtual ~efficiencyFunction();

        //- Disallow default bitwise assignment
        void operator=(const efficiencyFunction&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
