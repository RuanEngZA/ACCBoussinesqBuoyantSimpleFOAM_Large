{
    Istream& is4 = mesh.solutionDict().subDict("General").lookup("HeatOnOff");
    is4.format(IOstream::ASCII); 
    is4 >> heatOnOff;
  
  
  
    alphat = turbulence->nut()/Prt;
    alphat.correctBoundaryConditions();

    volScalarField alphaEff("alphaEff", turbulence->nu()/Pr + alphat);

    fvScalarMatrix TEqn
    (
        fvm::div(phi, T)
      - fvm::laplacian(alphaEff, T) - heatOnOff*energySource/(rho*cpa)
     ==
        fvOptions(T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    fvOptions.correct(T);

    rhok = 1.0 - beta*(T - TRef);
}
