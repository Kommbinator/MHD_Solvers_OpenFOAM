fvMesh& mesh = fluidRegions[i];
volVectorField& U = UFluid[i];
surfaceScalarField& phi = phiFluid[i];
volScalarField& p = pFluid[i];
volScalarField& PotE = PotEFluid[i];
volVectorField& B0 = B0Fluid[i];
volVectorField& lorentz = lorentzFluid[i];
volVectorField& J = JFluid[i];
volScalarField& divJ = divJFluid[i];
volScalarField& sigma = sigmaFluid[i];

//dimensionedScalar& sigma = sigmaFluid[i];
dimensionedScalar& rho = rhoFluid[i];
//dimensionedScalar& nu = nuFluid[i];

IOMRFZoneList& MRF = MRFfluid[i];
fv::options& fvOptions = fluidFvOptions[i];
    
const label pRefCell = pRefCellFluid[i];
const scalar pRefValue = pRefValueFluid[i];

const label PotERefCell = PotERefCellFluid[i];
const scalar PotERefValue = PotERefValueFluid[i];

incompressible::turbulenceModel& turbulence = turbulenceFluid[i];
