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

//Electric potential reference point
label PotERefCell = 0;
scalar PotERefValue = 0.0;
setRefCell(PotE, mesh.solutionDict().subDict("PIMPLE"), PotERefCell, PotERefValue);
mesh.setFluxRequired(PotE.name());
