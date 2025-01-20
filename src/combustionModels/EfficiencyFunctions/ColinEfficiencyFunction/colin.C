#include "colin.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcCurl.H"
#include "fvcLaplacian.H"
#include "DTF.H"

namespace Foam
{
namespace efficiencyFunctionModels
{
    defineTypeNameAndDebug(colin, 0);

    addToRunTimeSelectionTable
    (
        efficiencyFunction,
        colin,
        dictionary
    );

} // namespace efficiencyFunctionModels
} // namespace Foam

// Constructors
Foam::efficiencyFunctionModels::colin::colin
(
    const dictionary& dict,
    const fvMesh& mesh,
    const compressibleMomentumTransportModel& turb,
    scalar Fmax,
    const Foam::combustionModels::DTF& dtfModel
)
:
    efficiencyFunction(dict, mesh, turb, Fmax, dtfModel),
    coeffsDict_(dict.subDict("colinCoeffs")),
    alpha_(coeffsDict_.lookup<scalar>("alpha")),
    dtfModel_(dtfModel),
    uPrime_
    (
        IOobject
        (
            "uPrime",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimLength/dimTime, 0.0)
    ),
    sFilter_(mesh_)
{}


void Foam::efficiencyFunctionModels::colin::correct() {

    scalar delta_L = dtfModel_.deltaL();
    scalar SL_ = dtfModel_.SL();
    scalar Filters_ = dtfModel_.filters();
    const volScalarField& delta_ = mesh_.lookupObject<volScalarField>("delta");
    const volScalarField& Fc_ = dtfModel_.F();



    // Check validity of inputs
    if (delta_L <= 0) {
        FatalErrorInFunction << "Invalid flame thickness (delta_L): " << delta_L << exit(FatalError);
    }
    if (SL_ <= 0) {
        FatalErrorInFunction << "Invalid flame speed (SL_): " << SL_ << exit(FatalError);
    }

    Info << "Flame Thickness (delta_L): " << delta_L << endl;
    Info << "Flame Speed (SL_): " << SL_ << endl;

    // Calculate uPrime_ using filtered U field
    volVectorField U_hat(turb_.U());
    Info << "Uhat_ computed successfully." << endl;

    for (int i = 0; i < Filters_; i++) {
        U_hat = sFilter_(U_hat);
    }

    // Ensure delta_ is valid
    if (delta_.internalField().size() == 0) {
        FatalErrorInFunction << "delta_ field is empty or uninitialized." << exit(FatalError);
    }

    uPrime_ = 2.0 * mag(pow(delta_, 3) * fvc::laplacian(fvc::curl(U_hat)));



    

    // Calculate the efficiency
    forAll(E_, celli) {
        scalar delta_e = Fc_[celli] * delta_L;

        if (delta_e <= 0) {
            WarningInFunction << "Invalid delta_e at cell " << celli << endl;
            E_[celli] = 1.0; // Assign minimum efficiency to avoid invalid calculations
            continue;
        }

        scalar delta_r = delta_e / delta_L;
        scalar delta_R = delta_e / (delta_L * Fc_[celli]);
        
        if (delta_r <= 0 || delta_R <= 0) {
            FatalErrorInFunction << "Invalid delta_r or delta_R values." << exit(FatalError);
        }       

        scalar up_SL_r = (uPrime_[celli] * delta_e) / SL_;
        
        if (up_SL_r <= 0) {
            WarningInFunction << "Invalid up_SL_r at cell " << celli << endl;
            E_[celli] = 1.0; // Assign minimum efficiency to avoid invalid calculations
            continue;
        }

        scalar numerator = 1.0 + alpha_ * 0.75 * exp(-1.2 / pow(up_SL_r, 0.3)) * pow(delta_r, 2.0 / 3.0) * up_SL_r;
        scalar denominator = 1.0 + alpha_ * 0.75 * exp(-1.2 / pow(up_SL_r, 0.3)) * pow(delta_R, 2.0 / 3.0) * up_SL_r;

        if (denominator <= 0) {
            WarningInFunction << "Invalid denominator at cell " << celli << endl;
            E_[celli] = 1.0; // Assign minimum efficiency to avoid division by zero
            continue;
        }

        scalar efficiency = numerator / denominator;

        // Enforce the condition E >= 1
        E_[celli] = max(efficiency, 1.0);
    }
}

// Destructor
Foam::efficiencyFunctionModels::colin::~colin()
{}

