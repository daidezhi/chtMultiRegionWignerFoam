{
    Info<< "\nSolving for solid region "<< solidRegions[i].name() << endl;

    wignerEnergy.update();

    fvScalarMatrix hEqn
    (
        fvm::ddt(betav*rho, h)
      - (
            thermo.isotropic()
            ? fvm::laplacian(betav*thermo.alpha(), h, "laplacian(alpha,h)")
            : fvm::laplacian(betav*taniAlpha(), h, "laplacian(alpha,h)")
        )
      ==
        fvOptions(rho, h)
      + fvc::Sp(rho, wignerEnergy.SDot())   // Wigner energy source
    );

    hEqn.relax();

    fvOptions.constrain(hEqn);


    if (finalIter)
    {
        mesh.data().setFinalIteration(true);
    }

    hEqn.solve(h.select(finalIter));

    fvOptions.correct(h);

    thermo.correct();

    Info<< "Min/max T:" << min(thermo.T()).value() << ' '
        << max(thermo.T()).value() << endl;

    if (finalIter)
    {
        mesh.data().setFinalIteration(false);
    }
}
