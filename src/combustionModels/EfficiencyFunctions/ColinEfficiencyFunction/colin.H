#ifndef colin_H
#define colin_H

#include "efficiencyFunction.H"
#include "simpleFilter.H"
#include "momentumTransportModel.H"
 
namespace Foam
{
namespace efficiencyFunctionModels
{

/*---------------------------------------------------------------------------*\
                        Class colin Declaration
\*---------------------------------------------------------------------------*/

class colin
:
    public efficiencyFunction
{
    // Private Data
    dictionary coeffsDict_;
    scalar alpha_;
    const Foam::combustionModels::DTF& dtfModel_;
    volScalarField uPrime_;
   simpleFilter sFilter_;
    

public:

    //- Runtime type information
    TypeName("colin");

    // Constructors
    colin
    (
        const dictionary&,
        const fvMesh&,
        const compressibleMomentumTransportModel& turb,
        scalar F,
        const Foam::combustionModels::DTF& dtfModel
    );

    // Disallow default bitwise copy construction
    colin(const colin&) = delete;

    // Destructor
    virtual ~colin();

/*
    //- Access function to filter width
    inline const volScalarField& delta() const
        {
            return delta();
        }    
*/ 

    // Correct function
    virtual void correct();

    // Disallow default bitwise assignment
    void operator=(const colin&) = delete;
};

} // End namespace efficiencyFunctionModels
} // End namespace Foam

#endif
