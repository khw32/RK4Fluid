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

Application
    fsiFoam

Description
    Finite volume fluid structure interaction solver based on partitioned
    approach and strong coupling.

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "RK4fluidSolidInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicFvMesh.H"
#   include "createSolidMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    RK4fluidSolidInterface fsi(mesh, solidMesh);

    Info<< "\nStarting time loop\n" << endl;

    for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time = " << runTime.value()
            << " (dt = " << runTime.deltaT().value() << ")" << nl << endl;
        
        // Make the initialization by calling member function

        Info<< "\ninitializeFields\n" << endl;
        fsi.initializeFields();
         
Info<< "\nupdateInterpolator\n" << endl;
        fsi.updateInterpolator();

        scalar residualNorm = 0;

        // Make predictions of initial forces and evoke solid solver and update
        // the residual to prepare for the calculation loop

        if (fsi.predictor())
        {
Info<< "\nupdateForce\n" << endl;
            fsi.updateForce();
Info<< "\nsolid solver\n" << endl;
            fsi.stress().evolve();

            residualNorm =
                fsi.updateResidual();
        }

        // Loop the fluid-solid interaction procedure which is firstly updating
        // the solid displacement based on the solid solver and then transfer
        // the deformed mesh to fluid side followed by asking fluid solver to
        // calculate the pressure and velocity. Again the pressure and viscous
        // forces are updated and transported to the solid side to solve new
        // velocity and displacement

        do
        {
Info<< "\nIncrease outer corrector\n" << endl;
            fsi.outerCorr()++;

Info<< "\nupdate displacement\n" << endl;
            fsi.updateDisplacement();

Info<< "\nmove solid mesh\n" << endl;
            fsi.moveFluidMesh();
Info<< "\ncheck if statement\n" << endl;
            if (fsi.outerCorr() == 1)
            {
Info<< "\nCase_1 outer corrector = 1\n" << endl;
                fsi.flow().evolve();
            }                

            else

            {
Info<< "\nCase_2 outer corrector != 1\n" << endl;
                fsi.flow().evolveUpdate();
            } 

Info<< "\nUpdate Force\n" << endl;
            fsi.updateForce();
Info<< "\nSolve solid\n" << endl;
            fsi.stress().evolve();

Info<< "\nUpdate residual\n" << endl;
            residualNorm =
                fsi.updateResidual();
Info<< "\ncheck while statement\n" << endl;
        }

       // Check residual and iterate the loop until the requirement is met
        while
        (
            (residualNorm > fsi.outerCorrTolerance())
         && (fsi.outerCorr() < fsi.nOuterCorr())
        );

        // The final flow field and stress field are written to the output

Info<< "\nupdate flow fields\n" << endl;
        fsi.flow().updateFields();
Info<< "\nupdate total fields\n" << endl;
        fsi.stress().updateTotalFields();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //