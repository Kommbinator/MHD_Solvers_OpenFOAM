if (finalIter)
{
    mesh.data::add("finalIteration", true);
}

//Poisson equation for electric potential
fvScalarMatrix PotEEqn
(
fvm::laplacian(sigmaEMfield,PotE)
);

//Solving Poisson equation
PotEEqn.solve();
PotE.correctBoundaryConditions();

//Computation of current density at cell faces
surfaceScalarField jnSolid = -fvc::interpolate(sigmaEMfield) * fvc::snGrad(PotE)*mesh.magSf();	 

//Current density at face center for J visualisation or for the scheme (63) NI paper, part 2
surfaceVectorField jnvSolid = jnSolid * mesh.Cf();

//Interpolation of current density at cell center
J = fvc::surfaceIntegrate(jnvSolid) - (fvc::surfaceIntegrate(jnSolid) * mesh.C());	

//Update current density distribution and boundary condition
J.correctBoundaryConditions();


if (finalIter)
{
    mesh.data::remove("finalIteration");
}
