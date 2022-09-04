Info<< "Reading MHD transportProperties\n" << endl;
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
