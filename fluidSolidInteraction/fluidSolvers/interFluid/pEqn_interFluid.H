        volScalarField rAU_ = 1.0/UEqn.A();
        surfaceScalarField rAUf_ = fvc::interpolate(rAU_);

        U_ = rAU_*UEqn.H();
        surfaceScalarField phiU_("phiU", (fvc::interpolate(U_) & mesh.Sf()));
   
	if (pd_.needReference())
 	   {
               Info<< "pd needReference" << nl << endl;
 	       adjustPhi(phi_, U_, pd_);
 	   }

        phi_ = phiU_ +
    		(
        	fvc::interpolate(interface_.sigmaK())*fvc::snGrad(alpha1_)
      		- ghf_*fvc::snGrad(rho_)
    		)*rAUf_*mesh.magSf();


        //    while (pimple.correctNonOrthogonal())
        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
        	fvScalarMatrix pdEqn
        	(
        	    fvm::laplacian(rAUf_, pd_) == fvc::div(phi_)
        	);

        	pdEqn.setReference(pdRefCell, pdRefValue);

        	pdEqn.solve
        	(
          //  mesh.solutionDict().solver(pd.select(pimple.finalInnerIter())) pFinal ???
          // could be found at pimpleControlI.H    
        	);

       // if (pimple.finalNonOrthogonalIter())
        if (nonOrth == nNonOrthCorr)
       		 {
       		     phi_ -= pdEqn.flux();
       		 }
    	}

        U_ += rAU_*fvc::reconstruct((phi_ - phiU_)/rAUf_);
        U_.correctBoundaryConditions();
        gradU_ = fvc::grad(U_);

# include "continuityErrs_interFluid.H"

