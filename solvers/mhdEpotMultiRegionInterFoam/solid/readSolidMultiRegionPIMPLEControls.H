const dictionary& pimple = mesh.solutionDict().subDict("PIMPLE");

//Since the electrical potential equation is steady in the solid domain, the Pimple control tools are not used
//pimple.getOrDefault<int>("nNonOrthogonalCorrectors", 0);
