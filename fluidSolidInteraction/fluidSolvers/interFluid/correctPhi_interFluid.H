
# include "continuityErrs_interFluid.H"

    wordList pcorrTypes
    (
        pd_.boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    for (label i = 0; i < pd_.boundaryField().size(); i++)
    {
        if (pd_.boundaryField()[i].fixesValue())
        {
            pcorrTypes[i] = fixedValueFvPatchScalarField::typeName;
        }
    }


       volScalarField pcorr_
    (
        IOobject
        (
            "pcorr",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("pcorr", pd_.dimensions(), 0.0),
        pcorrTypes
    );

    dimensionedScalar rAUf_
    (
        "(1|A(U))",
        dimTime/rho_.dimensions(),
        runTime().deltaT().value()
    );

    phi_ = (fvc::interpolate(U_) & mesh.Sf());

    adjustPhi(phi_, U_, pcorr_); // adjust the balance of fluxes to obey continuity for cases 
    // which do not have a pressure boundary

    mesh.schemesDict().setFluxRequired(pcorr_.name());

    //while (pimple.correctNonOrthogonal())
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        // Solve for pcorr such that the divergence of the corrected flux
        // matches the divU provided (from previous iteration, time-step...)
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rAUf_, pcorr_) == fvc::div(phi_)
        );

        pcorrEqn.setReference(pdRefCell, pdRefValue);
        pcorrEqn.solve();

       // if (pimple.finalNonOrthogonalIter())
       if (nonOrth == nNonOrthCorr)
        {
            phi_ -= pcorrEqn.flux();
        }
    }
# include "continuityErrs_interFluid.H"
# include "CourantNo_interFluid.H"

    // recalculate rhoPhi from pho and corrected phi
    rhoPhi_ = fvc::interpolate(rho_)*phi_;

