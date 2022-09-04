// Initialise fluid field pointer lists
PtrList<volVectorField> UFluid(fluidRegions.size());
PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
PtrList<volScalarField> pFluid(fluidRegions.size());
PtrList<volScalarField> PotEFluid(fluidRegions.size());
PtrList<volVectorField> B0Fluid(fluidRegions.size());
PtrList<volVectorField> lorentzFluid(fluidRegions.size());
PtrList<volVectorField> JFluid(fluidRegions.size());
PtrList<volScalarField> divJFluid(fluidRegions.size());
PtrList<volScalarField> sigmaFluid(fluidRegions.size());


PtrList<IOdictionary> transportPropertiesFluid(fluidRegions.size());
PtrList<dimensionedScalar> rhoFluid(fluidRegions.size());
//PtrList<dimensionedScalar> sigmaFluid(fluidRegions.size());

PtrList<singlePhaseTransportModel> laminarTransportFluid(fluidRegions.size());
PtrList<incompressible::turbulenceModel> turbulenceFluid(fluidRegions.size());

List<label> pRefCellFluid(fluidRegions.size());
List<scalar> pRefValueFluid(fluidRegions.size());

List<label> PotERefCellFluid(fluidRegions.size());
List<scalar> PotERefValueFluid(fluidRegions.size());

PtrList<IOMRFZoneList> MRFfluid(fluidRegions.size());
PtrList<fv::options> fluidFvOptions(fluidRegions.size());

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
    Info<< "    Adding to pFluid\n" << endl;    
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
    
    sigmaFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "sigma",
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
 
    laminarTransportFluid.set(i, new singlePhaseTransportModel(UFluid[i], phiFluid[i]));    
    
    
    Info<< "    Adding to turbulenceFluid\n" << endl;
    turbulenceFluid.set
    (
        i,
        incompressible::turbulenceModel::New
        (
            UFluid[i],
            phiFluid[i],
            laminarTransportFluid[i]
        ).ptr()
    );
    Info<< "    Adding to transportPropertiesFluid\n" << endl;
    transportPropertiesFluid.set
    (
        i,        
        new IOdictionary 
        (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            fluidRegions[i],
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
        )
    );
    
    //Density
    rhoFluid.set
    (
        i,
        new dimensionedScalar
        (
            "rho",
            dimDensity,
            transportPropertiesFluid[i]
        )
    );
    
   /* sigmaFluid.set
    (
        i,
        new dimensionedScalar
        (
            "sigma",
            dimensionSet(-1, -3, 3, 0, 0, 2, 0),
            transportPropertiesFluid[i]
        )
    ); */   

    //Kinematic viscosity
   /* nuFluid.set
    (
        i,        
        new dimensionedScalar
        (
            "nu",
            dimViscosity,
            transportPropertiesFluid[i]
        )
    );*/
    

    pRefCellFluid[i] = 0;
    pRefValueFluid[i] = 0.0;
    
    setRefCell(pFluid[i], fluidRegions[i].solutionDict().subDict("PIMPLE"), pRefCellFluid[i], pRefValueFluid[i]);
    fluidRegions[i].setFluxRequired(pFluid[i].name());
    
    PotERefCellFluid[i] = 0;
    PotERefValueFluid[i] = 0.0;
    
    setRefCell(PotEFluid[i], fluidRegions[i].solutionDict().subDict("PIMPLE"), PotERefCellFluid[i], PotERefValueFluid[i]);
    fluidRegions[i].setFluxRequired(PotEFluid[i].name());
    
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
