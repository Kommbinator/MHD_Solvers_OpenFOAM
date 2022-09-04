// Initialise fluid field pointer lists
PtrList<volVectorField> UFluid(fluidRegions.size());
PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
PtrList<volScalarField> p_rghFluid(fluidRegions.size());
PtrList<volScalarField> PotEFluid(fluidRegions.size());
PtrList<volVectorField> B0Fluid(fluidRegions.size());
PtrList<volVectorField> lorentzFluid(fluidRegions.size());
PtrList<volVectorField> JFluid(fluidRegions.size());
PtrList<volScalarField> divJFluid(fluidRegions.size());

PtrList<volScalarField> rhoFluid(fluidRegions.size());
PtrList<surfaceScalarField> rhoPhiFluid(fluidRegions.size());
PtrList<volScalarField> sigmaEMfieldFluid(fluidRegions.size());
PtrList<uniformDimensionedScalarField> hRefFluid(fluidRegions.size());
PtrList<volScalarField> ghFluid(fluidRegions.size());
PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());
PtrList<volScalarField> pFluid(fluidRegions.size());

PtrList<immiscibleIncompressibleTwoPhaseMhdMixture> mixtureFluid(fluidRegions.size());
PtrList<incompressible::turbulenceModel> turbulenceFluid(fluidRegions.size());

List<label> pRefCellFluid(fluidRegions.size());
List<scalar> pRefValueFluid(fluidRegions.size());

List<label> PotERefCellFluid(fluidRegions.size());
List<scalar> PotERefValueFluid(fluidRegions.size());

PtrList<surfaceScalarField> alphaPhiUnFluid(fluidRegions.size());
PtrList<const bool> alphaRestartFluid(fluidRegions.size());
PtrList<surfaceScalarField> alphaPhi10Fluid(fluidRegions.size());


PtrList<IOMRFZoneList> MRFfluid(fluidRegions.size());
PtrList<fv::options> fluidFvOptions(fluidRegions.size());

Info<< "\nReading g" << endl;
const uniformDimensionedVectorField& g = meshObjects::gravity::New(runTime);

// Populate fluid field pointer lists
forAll(fluidRegions, i)
{   
    
    Info<< "    Adding to UFluid\n" << endl;
    UFluid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );  
    Info<< "    Adding to p_rghFluid\n" << endl;    
    p_rghFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p_rgh",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );    
    
    Info<< "    Adding to phiFluid\n" << endl;
    phiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            linearInterpolate(UFluid[i])
                & fluidRegions[i].Sf()
        )
    );    
    
    Info<< "    Adding to PotEFluid\n" << endl;    
    PotEFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "PotE",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );   
    
    Info<< "    Adding to B0Fluid\n" << endl;
    B0Fluid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "B0",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );    
    
    lorentzFluid.set
    (
        i,
        new volVectorField
        (
        IOobject
        (
            "lorentz",
            runTime.timeName(),
            fluidRegions[i],
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluidRegions[i],
        dimensionedVector("lorentz", dimensionSet(1,-2,-2,0,0,0,0), Foam::vector(0,0,0))
        )
    );    
        
    JFluid.set
    (
        i,
        new volVectorField
        (
        IOobject
        (
            "J",
            runTime.timeName(),
            fluidRegions[i],
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluidRegions[i],
        dimensionedVector("J", dimensionSet(0,-2,0,0,0,1,0), Foam::vector(0,0,0))
        )
    );

        //Divergence of J
    divJFluid.set
    (
        i,
        new volScalarField
        (
        IOobject
        (
            "divJ",
            runTime.timeName(),
            fluidRegions[i],
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluidRegions[i],
        dimensionedScalar("divJ", dimensionSet(0,-3,0,0,0,1,0), Foam::scalar(0))
        )
    );       
 
        
    Info<< "Reading transportProperties\n" << endl;
    
    mixtureFluid.set(i, new immiscibleIncompressibleTwoPhaseMhdMixture(UFluid[i], phiFluid[i])); 
    
    volScalarField& alpha1(mixtureFluid[i].alpha1());
    volScalarField& alpha2(mixtureFluid[i].alpha2());

    const dimensionedScalar& rho1 = mixtureFluid[i].rho1();
    const dimensionedScalar& rho2 = mixtureFluid[i].rho2();

    const dimensionedScalar& sigma1 = mixtureFluid[i].sigma1();
    const dimensionedScalar& sigma2 = mixtureFluid[i].sigma2();    
    
    Info<< "   Adding rhoFluid\n" << endl;    
    
    rhoFluid.set
    (   i,
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            alpha1*rho1 + alpha2*rho2
        )
    ); 
    
    rhoPhiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "rhoPhi",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::interpolate(rhoFluid[i])*phiFluid[i]
        )
    );   
    
    sigmaEMfieldFluid.set
    (   i,
        new volScalarField
        (
            IOobject
            (
                "sigmaEMfield",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            alpha1*sigma1 + alpha2*sigma2
        )
    );    
    
    
Info<< "    Adding to hRefFluid\n" << endl;
    hRefFluid.set
    (
        i,
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "hRef",
                runTime.constant(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedScalar("hRef", dimLength, Zero) // uses name
        )
    );

    dimensionedScalar ghRef
    (
        mag(g.value()) > SMALL
      ? g & (cmptMag(g.value())/mag(g.value()))*hRefFluid[i]
      : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
    );

    Info<< "    Adding to ghFluid\n" << endl;
    ghFluid.set
    (
        i,
        new volScalarField
        (
            "gh",
            (g & fluidRegions[i].C()) - ghRef
        )
    );

    Info<< "    Adding to ghfFluid\n" << endl;
    ghfFluid.set
    (
        i,
        new surfaceScalarField
        (
            "ghf",
            (g & fluidRegions[i].Cf()) - ghRef
        )
    );    
    
    pFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            p_rghFluid[i] + rhoFluid[i]*ghFluid[i]
        )
    );    
           
    
    Info<< "    Adding to turbulenceFluid\n" << endl;
    turbulenceFluid.set
    (
        i,
        incompressible::turbulenceModel::New
        (
            UFluid[i],
            phiFluid[i],
            mixtureFluid[i]
        ).ptr()
    );    

    pRefCellFluid[i] = 0;
    pRefValueFluid[i] = 0.0;
    
    setRefCell(pFluid[i], p_rghFluid[i], fluidRegions[i].solutionDict().subDict("PIMPLE"), pRefCellFluid[i], pRefValueFluid[i]);
    fluidRegions[i].setFluxRequired(p_rghFluid[i].name());
    
    if (p_rghFluid[i].needReference())
    {
        pFluid[i] += dimensionedScalar
        (
            "pFluid",
            pFluid[i].dimensions(),
            pRefValueFluid[i] - getRefCellValue(pFluid[i], pRefCellFluid[i])
        );
        p_rghFluid[i] = pFluid[i] - rhoFluid[i]*ghFluid[i];
    }    
    
    fluidRegions[i].setFluxRequired(alpha1.name());

    
    PotERefCellFluid[i] = 0;
    PotERefValueFluid[i] = 0.0;
    
    setRefCell(PotEFluid[i], fluidRegions[i].solutionDict().subDict("PIMPLE"), PotERefCellFluid[i], PotERefValueFluid[i]);
    fluidRegions[i].setFluxRequired(PotEFluid[i].name());
    
    // MULES compressed flux is registered in case scalarTransport FO needs it.
    alphaPhiUnFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "alphaPhiUn",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
        fluidRegions[i],
        dimensionedScalar(phiFluid[i].dimensions(), Zero)
        )
    );  
    
    //Alpha fluxes
    IOobject alphaPhi10Header
    (
        IOobject::groupName("alphaPhi0", alpha1.group()),
        runTime.timeName(),
        fluidRegions[i],
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    );

    alphaRestartFluid.set
    (
        i,
        new const bool 
        (
            alphaPhi10Header.typeHeaderOk<surfaceScalarField>(true)
        )
    );

    if (alphaRestartFluid[i])
    {
        Info << "   Restarting alpha!" << endl;
    }

    // MULES flux from previous time-step
    alphaPhi10Fluid.set
    (
        i, 
        new surfaceScalarField
        (
            alphaPhi10Header,
            phiFluid[i]*fvc::interpolate(alpha1)
        )
    );    
    
    Info<< "    Adding MRF\n" << endl;
    MRFfluid.set
    (
        i,
        new IOMRFZoneList(fluidRegions[i])
    );

    Info<< "    Adding fvOptions\n" << endl;
    fluidFvOptions.set
    (
        i,
        new fv::options(fluidRegions[i])
    );   
    
}
