fvMesh& mesh = fluidRegions[i];

immiscibleIncompressibleTwoPhaseMhdMixture& mixture  = mixtureFluid[i];

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());
const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();
const dimensionedScalar& sigma1 = mixture.sigma1();
const dimensionedScalar& sigma2 = mixture.sigma2();

volVectorField& U = UFluid[i];
surfaceScalarField& phi = phiFluid[i];
volScalarField& p = pFluid[i];
volScalarField& p_rgh = p_rghFluid[i];
volScalarField& PotE = PotEFluid[i];
volVectorField& B0 = B0Fluid[i];
volVectorField& lorentz = lorentzFluid[i];
volVectorField& J = JFluid[i];
volScalarField& divJ = divJFluid[i];

volScalarField& rho = rhoFluid[i];
surfaceScalarField& rhoPhi = rhoPhiFluid[i];
volScalarField& sigmaEMfield = sigmaEMfieldFluid[i];

const volScalarField& gh = ghFluid[i];
const surfaceScalarField& ghf = ghfFluid[i];

const bool alphaRestart = alphaRestartFluid[i];
surfaceScalarField& alphaPhi10 = alphaPhi10Fluid[i];

// MULES Correction
tmp<surfaceScalarField> talphaPhi1Corr0;
    
IOMRFZoneList& MRF = MRFfluid[i];
fv::options& fvOptions = fluidFvOptions[i];
    
const label pRefCell = pRefCellFluid[i];
const scalar pRefValue = pRefValueFluid[i];

const label PotERefCell = PotERefCellFluid[i];
const scalar PotERefValue = PotERefValueFluid[i];

incompressible::turbulenceModel& turbulence = turbulenceFluid[i];
