/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "RK4Fluid.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "findRefCell.H"
#include "adjustPhi.H"
#include "fluidSolidInterface.H"
#include "fixedGradientFvPatchFields.H"
#include "fixedFluxPressureFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace RK4fluidSolvers
    {
        // * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

        defineTypeNameAndDebug(RK4Fluid, 0);
        addToRunTimeSelectionTable(RK4fluidSolver, RK4Fluid, dictionary);

        // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

        RK4Fluid::RK4Fluid(const fvMesh& mesh)
        :
        RK4fluidSolver(this->typeName, mesh),
        U_
        (
            IOobject
            (
                "U",
                runTime().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        ),
        Uold_
        (
            IOobject
            (
                "Uold",
                runTime().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            U_
        ),
        Uc_
        (
            IOobject
            (
                "Uc",
                runTime().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            U_
        ),
        dU_
        (
            IOobject
            (
                "dU",
                runTime().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            U_
        ),
        p_
        (
            IOobject
            (
                "p",
                runTime().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        ),
        gradp_(fvc::grad(p_)),
        gradU_(fvc::grad(U_)),
        phi_
        (
            IOobject
            (
                "phi",
                runTime().timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::interpolate(U_) & mesh.Sf()
        ),
        laminarTransport_(U_, phi_),
        turbulence_
        (
            incompressible::turbulenceModel::New
            (
                U_, phi_, laminarTransport_
            )
        ),
        rho_
        (
            IOdictionary
            (
                IOobject
                (
                    "transportProperties",
                    runTime().constant(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            ).lookup("rho")
        )
    {}

        // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

        const volVectorField& RK4Fluid::U() const
        {
            return U_;
        }

        const volVectorField& RK4Fluid::Uold() const
        {
            return Uold_;
        }

        const volVectorField& RK4Fluid::Uc() const
        {
            return Uc_;
        }

        const volVectorField& RK4Fluid::dU() const
        {
            return dU_;
        }

        const volScalarField& RK4Fluid::p() const
        {
            return p_;
        }


        //- Patch viscous force (N/m2)
        tmp<vectorField> RK4Fluid::patchViscousForce(const label patchID) const
        {
            tmp<vectorField> tvF
            (
                new vectorField(mesh().boundary()[patchID].size(), vector::zero)
            );

            tvF() =
                rho_.value()
            *(
                    mesh().boundary()[patchID].nf()
                & turbulence_->devReff()().boundaryField()[patchID]
                );

            //     vectorField n = mesh().boundary()[patchID].nf();
            //     tvF() -= n*(n&tvF());

            return tvF;
        }

        //- Patch pressure force (N/m2)
        tmp<scalarField> RK4Fluid::patchPressureForce(const label patchID) const
        {
            tmp<scalarField> tpF
            (
                new scalarField(mesh().boundary()[patchID].size(), 0)
            );

            tpF() = rho_.value()*p().boundaryField()[patchID];

            return tpF;
        }

        //- Patch viscous force (N/m2)
        tmp<vectorField> RK4Fluid::faceZoneViscousForce
        (
            const label zoneID,
            const label patchID
        ) const
        {
            vectorField pVF = patchViscousForce(patchID);

            tmp<vectorField> tvF
            (
                new vectorField(mesh().faceZones()[zoneID].size(), vector::zero)
            );
            vectorField& vF = tvF();

            const label patchStart =
                mesh().boundaryMesh()[patchID].start();

            forAll(pVF, i)
            {
                vF[mesh().faceZones()[zoneID].whichFace(patchStart + i)] =
                    pVF[i];
            }

            // Parallel data exchange: collect pressure field on all processors
            reduce(vF, sumOp<vectorField>());


            return tvF;
        }

        //- Patch pressure force (N/m2)
        tmp<scalarField> RK4Fluid::faceZonePressureForce
        (
            const label zoneID,
            const label patchID
        ) const
        {
            scalarField pPF = patchPressureForce(patchID);

            tmp<scalarField> tpF
            (
                new scalarField(mesh().faceZones()[zoneID].size(), 0)
            );
            scalarField& pF = tpF();

            const label patchStart =
                mesh().boundaryMesh()[patchID].start();

            forAll(pPF, i)
            {
                pF[mesh().faceZones()[zoneID].whichFace(patchStart + i)] =
                    pPF[i];
            }

            // Parallel data exchange: collect pressure field on all processors
            reduce(pF, sumOp<scalarField>());

            return tpF;
        }

        tmp<scalarField> RK4Fluid::faceZoneMuEff
        (
            const label zoneID,
            const label patchID
        ) const
        {
            scalarField pMuEff =
                rho_.value()*turbulence_->nuEff()().boundaryField()[patchID];

            tmp<scalarField> tMuEff
            (
                new scalarField(mesh().faceZones()[zoneID].size(), 0)
            );
            scalarField& muEff = tMuEff();

            const label patchStart =
                mesh().boundaryMesh()[patchID].start();

            forAll(pMuEff, i)
            {
                muEff[mesh().faceZones()[zoneID].whichFace(patchStart + i)] =
                    pMuEff[i];
            }

            // Parallel data exchange: collect pressure field on all processors
            reduce(muEff, sumOp<scalarField>());

            return tMuEff;
        }

        void RK4Fluid::evolve()
        {
            Info << "Evolving fluid solver" << endl;

            const fvMesh& mesh = RK4fluidSolver::mesh();
            /*
            int nCorr(readInt(fluidProperties().lookup("nCorrectors")));

            int nNonOrthCorr =
                readInt(fluidProperties().lookup("nNonOrthogonalCorrectors"));
            */

            // Rung-Kutta Coefficients

            scalarList beta(4);

            beta[0] = 0.166666666667;
            beta[1] = 0.333333333333;
            beta[2] = 0.333333333333;
            beta[3] = 0.166666666667;

            Info << "\n b1 = " <<beta[0]<< "\n b2 = " <<beta[1]<<"\n b3 = " <<beta[2]<<"\n b4 = " <<beta[3]<< endl;
            /*
            // Create time controls

            bool adjustTimeStep =
                runTime.controlDict().lookupOrDefault("adjustTimeStep", false);

            scalar maxCo =
                runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);

            scalar maxDeltaT =
                runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);
            */
            // Prepare for the pressure solution

            label pRefCell = 0;
            scalar pRefValue = 0.0;
            setRefCell(p_, fluidProperties(), pRefCell, pRefValue);

            // Create the Poisson Matrix

            Info << "Creating poisson matrix" << endl;
                
            fvScalarMatrix pEqn
            (
                fvm::laplacian(p_)
            );
            pEqn.setReference(pRefCell, pRefValue);

            //     if(mesh.moving())
            //     {
            //         // Make the fluxes relative
            //         phi_ -= fvc::meshPhi(U_);
            //     }

            // Make the fluxes relative to the mesh motion

            fvc::makeRelative(phi_, U_);
            /*
            // Read Time Controls

            adjustTimeStep =
                runTime.controlDict().lookupOrDefault("adjustTimeStep", false);

            maxCo =
                runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);

            maxDeltaT =
                runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);
            */
            // CourantNo
            {
                scalar CoNum = 0.0;
                scalar meanCoNum = 0.0;
                scalar velMag = 0.0;

                if (mesh.nInternalFaces())
                {
                    surfaceScalarField SfUfbyDelta =
                        mesh.surfaceInterpolation::deltaCoeffs()*mag(phi_);

                    CoNum = max(SfUfbyDelta/mesh.magSf())
                        .value()*runTime().deltaT().value();

                    meanCoNum = (sum(SfUfbyDelta)/sum(mesh.magSf()))
                        .value()*runTime().deltaT().value();

                    velMag = max(mag(phi_)/mesh.magSf()).value();
                }

                Info<< "Courant Number mean: " << meanCoNum
                    << " max: " << CoNum
                    << " velocity magnitude: " << velMag << endl;
            }

            // Adjust DeltaT (Might not be required)
            /*
            if (adjustTimeStep)
            {
                scalar maxDeltaTFact = maxCo/(CoNum + SMALL);
                scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

                runTime.setDeltaT
                (
                    min
                    (
                        deltaTFact*runTime.deltaT().value(),
                        maxDeltaT
                    )
                );

                Info<< "deltaT = " <<  runTime.deltaT().value() << endl;
            }
            */

            // 4th-order Rung-Kutta Algorithm

            // Prepare working variables

            dimensionedScalar normalizedTime("normalizedTime", dimensionSet(0,0,1,0,0,0,0),SMALL);

            Uold_ = U_;
            Uc_ = U_;
            normalizedTime = runTime().deltaT();

            forAll (beta, i)
            {

                // Calculate flux

                phi_ = (fvc::interpolate(U_) & mesh.Sf());

                // Calculate U increment

                dU_ = normalizedTime * (fvc::laplacian(turbulence_->nuEff(), U_) - fvc::div(phi_,U_));

                // Update variables

                Uc_ = Uc_ + beta[i] * dU_;
                U_ = Uold_ + 0.5 * dU_;

                // Continuity error

                {
                volScalarField contErr = fvc::div(phi_);

                scalar sumLocalContErr = runTime().deltaT().value()*
                    mag(contErr)().weightedAverage(mesh.V()).value();

                scalar globalContErr = runTime().deltaT().value()*
                    contErr.weightedAverage(mesh.V()).value();

                Info<< "time step continuity errors : sum local = "
                    << sumLocalContErr << ", global = " << globalContErr << endl;
                }

                // Pressure correction

                U_.correctBoundaryConditions();
                solve(pEqn == fvc::div(U_)/runTime().deltaT());
                U_ = U_ - (fvc::grad(p_) * runTime().deltaT());
                U_.correctBoundaryConditions();

            }
            
            phi_ = (fvc::interpolate(U_) & mesh.Sf());
            turbulence_->correct();

            // Make the fluxes absolut to the mesh motion
            fvc::makeAbsolute(phi_, U_);

        }

        void RK4Fluid::evolveUpdate()
        {
            Info << "Solving fluid subiteration loop" << endl;

            const fvMesh& mesh = RK4fluidSolver::mesh();

            // Prepare for the pressure solution
            label pRefCell = 0;
            scalar pRefValue = 0.0;
            setRefCell(p_, fluidProperties(), pRefCell, pRefValue);
            
            Info << "Creating poisson matrix" << endl;

            int nCorr(readInt(fluidProperties().lookup("nCorrectors")));

            for (int corr = 0; corr < nCorr; corr++)  
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(p_)
                );
                pEqn.setReference(pRefCell, pRefValue);

                phi_ = (fvc::interpolate(U_) & mesh.Sf());

                // Pressure correction

                U_.correctBoundaryConditions();
                solve(pEqn == fvc::div(U_)/runTime().deltaT());
                //solve(pEqn == fvc::div(U_));
                U_ = U_ - (fvc::grad(p_) * runTime().deltaT());
                //U_ = U_ - (fvc::grad(p_));
                // U_.correctBoundaryConditions();

                // phi_ = (fvc::interpolate(U_) & mesh.Sf());
                phi_ -= pEqn.flux()* runTime().deltaT();
            
            // Continuity error

                {
                volScalarField contErr = fvc::div(phi_);

                scalar sumLocalContErr = runTime().deltaT().value()*
                    mag(contErr)().weightedAverage(mesh.V()).value();

                scalar globalContErr = runTime().deltaT().value()*
                    contErr.weightedAverage(mesh.V()).value();

                Info<< "time step continuity errors : sum local = "
                    << sumLocalContErr << ", global = " << globalContErr << endl;
                }

                // Make the fluxes relative to the mesh motion
                fvc::makeRelative(phi_, U_);
                U_.correctBoundaryConditions();

            }

                turbulence_->correct();

                // Make the fluxes absolut to the mesh motion
                fvc::makeAbsolute(phi_, U_);

        }

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    } // End namespace RK4fluidSolvers
} // End namespace Foam

// ************************************************************************* //
