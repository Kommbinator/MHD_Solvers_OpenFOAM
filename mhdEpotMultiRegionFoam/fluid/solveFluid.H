if (finalIter)
{
    mesh.data::add("finalIteration", true);
}

#include "UEqn.H"

// --- PISO loop
for (int corr=0; corr<nCorr; corr++)
{
    #include "pEqn.H"
}

turbulence.correct();

#include "PotEEqn.H"


if (finalIter)
{
    mesh.data::remove("finalIteration");
}
