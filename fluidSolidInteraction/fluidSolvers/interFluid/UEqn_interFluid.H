   // #include Uqn.H
    surfaceScalarField muEff_
    (
        "muEff",
        twoPhaseProperties_.muf()
      + fvc::interpolate(rho_*turbulence_->nut())
    );

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho_, U_)
      + fvm::div(rhoPhi_, U_)
      - fvm::laplacian(muEff_, U_)
      - (fvc::grad(U_) & fvc::grad(muEff_))
    //- fvc::div(muEff*(fvc::interpolate(dev(fvc::grad(U))) & mesh.Sf()))
    );

    UEqn.relax();

    //    if (pimple.momentumPredictor()) // a member function in solutionControl return momentumPredictor_
    // pimpleControl is a subClass of solutionControl
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    fvc::interpolate(interface_.sigmaK())*fvc::snGrad(alpha1_)
                  - ghf_*fvc::snGrad(rho_)
                  - fvc::snGrad(pd_)
                )*mesh.magSf()
            )
        );
    }
