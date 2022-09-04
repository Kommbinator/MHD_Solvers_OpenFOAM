Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

//Density
dimensionedScalar rho
(
    "rho",
    dimDensity,
    transportProperties
);

//Kinematic viscosity
dimensionedScalar nu
(
    "nu",
    dimViscosity,
    transportProperties
);

//Electrical conductivity
dimensionedScalar sigma
(
    "sigma",
    dimensionSet(-1, -3, 3, 0, 0, 2, 0),
    transportProperties
);
 

//Pressure field
Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//Velocity field
Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

Info<< "Reading field PotE\n" << endl;
volScalarField PotE
(
    IOobject
    (
        "PotE",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//Magnetic field
Info<< "Reading field B0\n" << endl;
volVectorField B0
(
    IOobject
    (
        "B0",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

    //Lorentz field
    volVectorField lorentz
    (
        IOobject
        (
            "lorentz",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("lorentz", dimensionSet(1,-2,-2,0,0,0,0), Foam::vector(0,0,0))
    );
    
    //J field
    volVectorField J
    (
        IOobject
        (
            "J",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("J", dimensionSet(0,-2,0,0,0,1,0), Foam::vector(0,0,0))
    );
    
        volScalarField divJ
    (
        IOobject
        (
            "divJ",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("divJ", dimensionSet(0,-3,0,0,0,1,0), Foam::scalar(0))
    );

//Pressure reference point
label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

//Electric potential reference point
label PotERefCell = 0;
scalar PotERefValue = 0.0;
setRefCell(PotE, mesh.solutionDict().subDict("PIMPLE"), PotERefCell, PotERefValue);
mesh.setFluxRequired(PotE.name());

//Turbulent part
singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
#include "createFvOptions.H"
