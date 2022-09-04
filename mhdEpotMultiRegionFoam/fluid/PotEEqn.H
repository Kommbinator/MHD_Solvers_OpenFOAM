//Solve potential equation

//Interpolating cross product u x B over mesh faces
surfaceScalarField psiub = fvc::interpolate(sigma) * ((fvc::interpolate(U) ^ fvc::interpolate(B0)) & mesh.Sf());

//Poisson equation for electric potential
fvScalarMatrix PotEEqn
(
fvm::laplacian(sigma,PotE)
==
fvc::div(psiub)
);

//Reference potential
PotEEqn.setReference(PotERefCell, PotERefValue);
//PotEEqn.setReference(PotERefCell, getRefCellValue(PotE, PotERefCell));

//Solving Poisson equation
PotEEqn.solve();
PotE.correctBoundaryConditions();

//Computation of current density at cell faces
surfaceScalarField jn = -fvc::interpolate(sigma) * fvc::snGrad(PotE)*mesh.magSf()  + psiub;	 


//Current density at face center for J visualisation or for the scheme (63) NI paper, part 2
surfaceVectorField jnv = jn * mesh.Cf();

//Current density at face center
surfaceVectorField jnvb = jn * (fvc::interpolate(B0)^mesh.Cf());

//Interpolation of current density at cell center
J = fvc::surfaceIntegrate(jnv) - (fvc::surfaceIntegrate(jn) * mesh.C());	

//---Lorentz force computation, scheme (62) NI paper, part 2
lorentz = - fvc::surfaceIntegrate(jnvb) - (mesh.C()^(fvc::surfaceIntegrate(jn*(fvc::interpolate(B0)))));

//Update current density distribution and boundary condition
//J.correctBoundaryConditions();

//lorentz = (J ^ B0);
lorentz.correctBoundaryConditions();

//Calculate the error of div J
divJ = fvc::div(jn);
