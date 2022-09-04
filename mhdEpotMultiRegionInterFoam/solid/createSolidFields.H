    // Initialise solid field pointer lists
PtrList<volScalarField> PotESolid(solidRegions.size());
PtrList<volVectorField> JSolid(solidRegions.size());
PtrList<volScalarField> sigmaEMfieldSolid(solidRegions.size());

//PtrList<IOdictionary> transportPropertiesSolid(fluidRegions.size());
//PtrList<dimensionedScalar> sigmaSolid(fluidRegions.size());

// Populate solid field pointer lists
forAll(solidRegions, i)
{
    
Info<< "    Adding to PotESolid\n" << endl;    
PotESolid.set
(
        i,
        new volScalarField
        (
            IOobject
            (
                "PotE",
                runTime.timeName(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i]
        )
    );   

    JSolid.set
    (
        i,
        new volVectorField
        (
        IOobject
        (
            "J",
            runTime.timeName(),
            solidRegions[i],
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        solidRegions[i],
        dimensionedVector("J", dimensionSet(0,-2,0,0,0,1,0), Foam::vector(0,0,0))
        )
    );
    
    sigmaEMfieldSolid.set
    (   i,
        new volScalarField
        (
            IOobject
            (
                "sigmaEMfield",
                runTime.timeName(),
                solidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i]
        )
    );         

}
