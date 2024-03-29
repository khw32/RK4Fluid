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

Class
    RK4Fluid

Description
    icoFoam fluid solver

Author
    Zeljko Tukovic, FSB.  All rights reserved.

SourceFiles
    RK4Fluid.C

\*---------------------------------------------------------------------------*/

#ifndef RK4Fluid_H
#define RK4Fluid_H

#include "RK4fluidSolver.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RK4fluidSolvers
{

/*---------------------------------------------------------------------------*\
                        Class RK4Fluid Declaration
\*---------------------------------------------------------------------------*/

class RK4Fluid
:
    public RK4fluidSolver
{
    // Private data

        //- Velocity field
        volVectorField U_;

        //- Previous Iteration Variables

        volVectorField Uold_;
        volVectorField Uc_;
        volVectorField dU_;

        //- Pressure field
        volScalarField p_;

        //- Pressure field
        volVectorField gradp_;

        //- Velocity gradient
        volTensorField gradU_;

        //- Flux field
        surfaceScalarField phi_;

        //- Transport model
        singlePhaseTransportModel laminarTransport_;

        //- Turbulence model
        autoPtr<incompressible::turbulenceModel> turbulence_;

        //- Density
        dimensionedScalar rho_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        RK4Fluid(const RK4Fluid&);

        //- Disallow default bitwise assignment
        void operator=(const RK4Fluid&);


public:

    //- Runtime type information
    TypeName("RK4Fluid");

    // Constructors

        //- Construct from components
        RK4Fluid(const fvMesh& mesh);

//         //- Construct from components
//         RK4Fluid
//         (
//             const word& type,
//             const fvMesh& mesh
//         );

    // Destructor

        virtual ~RK4Fluid()
        {}


    // Member Functions

        // Access

            //- Return velocity field
            virtual const volVectorField& U() const;

            //- Return velocity fields
            volVectorField& U()
            {
                return U_;
            }

            virtual const volVectorField& Uold() const;

            volVectorField& Uold()
            {
                return Uold_;
            }

            virtual const volVectorField& Uc() const;

            volVectorField& Uc()
            {
                return Uc_;
            }

            virtual const volVectorField& dU() const;

            volVectorField& dU()
            {
                return dU_;
            }

            //- Return pressure field
            virtual const volScalarField& p() const;

            //- Return pressure field
            volScalarField& p()
            {
                return p_;
            }

            //- Return pressure gradient
            volVectorField& gradp()
            {
                return gradp_;
            }

            //- Return flux field
            surfaceScalarField& phi()
            {
                return phi_;
            }

//             //- Return kinematic viscosity
//             const dimensionedScalar& nu() const
//             {
//                 return nu_;
//             }

            //- Density
            const dimensionedScalar& rho()
            {
                return rho_;
            }

            //- Patch viscous force (N/m2)
            virtual tmp<vectorField> patchViscousForce
            (
                const label patchID
            ) const;

            //- Patch pressure force (N/m2)
            virtual tmp<scalarField> patchPressureForce
            (
                const label patchID
            ) const;

            //- Patch viscous force (N/m2)
            virtual tmp<vectorField> faceZoneViscousForce
            (
                const label zoneID,
                const label patchID
            ) const;

            //- Patch pressure force (N/m2)
            virtual tmp<scalarField> faceZonePressureForce
            (
                const label zoneID,
                const label patchID
            ) const;

            //- Face zone effective dynamic viscosity
            virtual tmp<scalarField> faceZoneMuEff
            (
                const label zoneID,
                const label patchID
            ) const;

        // Edit
        
            //- Evolve the fluid solver
            virtual void evolve();
            virtual void evolveUpdate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RK4fluidSolvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
