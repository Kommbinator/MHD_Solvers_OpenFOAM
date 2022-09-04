// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
        fvm::ddt(U) + fvm::div(phi, U)
    + MRF.DDt(U)
    + turbulence.divDevReff(U)
    - (1.0/rho) * lorentz //Lorentz force term  
    ==
        fvOptions(U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (momentumPredictor)
{
    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
}
