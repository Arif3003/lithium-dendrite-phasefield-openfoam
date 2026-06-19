/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    phaseFieldElectroMechano

Description
    Transient segregated finite-volume solver of linear-elastic,
    small-strain deformation of a solid body, with optional thermal
    diffusion and thermal stresses.

    Simple linear elasticity structural analysis code with included phaseFieldElectro.
    Solves for the displacement vector field D(also generating the
    stress tensor field sigma), phase-field variable (incorporating surface anisotropy), chemical potential and electric potential using a Grand Potential formulation

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
// #include "solidDisplacementThermo.H"
#include "Switch.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "Random.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    //#include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    simpleControl simple(mesh);
    #include "createControls.H"
    // #include "createControl.H"
    #include "createFields.H"
    // #include "createFieldRefs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating displacement field\n" << endl;
    

    volVectorField q=dimx*fvc::grad(eta);
    volScalarField magq = 0.0*eta;
    // volScalarField T = 0.0*phi + initial;
    // volScalarField sumLHS = 0.0*phi;
    volVectorField q_4 = 0.0*q;
    volVectorField q_6 = 0.0*q;
    volScalarField ac_01 = 0.0*eta;
    volVectorField dAdq01 = eta*vector(0,0,0);
    volVectorField dadgradPhi = q*0.0;  //<--- phi is eta

    while (runTime.loop())
    {
       Info<< "Iteration: " << runTime.value() << nl << endl;
       Info<< "Time = " << runTime.timeName() << nl << endl;

       #include "readSolidDisplacementFoamControls.H"

       int iCorr = 0;
       scalar initialResidual = 0; 
        // forAll(eta, cellI) 
        //               {
        //                 etaOld[cellI]=eta[cellI];

        //               }  
        while (simple.correctNonOrthogonal())
        {
            label curTimeIndex = mesh.time().timeIndex();
    for (int i=0; i<30; i++)
        {
		if(curTimeIndex == 0)
				{
                 //! The interpolation equation is used for smoothening of the
                 //! phase-field variable
                                    #include "preAllen.H"
                                }



	}

            //    mesh.update();
                //! Solving the phase-field and chemical potential equations after
                //! updating the mesh
                Info<< "\nCalculating eta,mu and phi distribution\n" << endl;
                #include "etaEqn.H"
    }

       do
       {
            // if (thermo.thermalStress())
            // {
            //     volScalarField& T = thermo.T();
            //     fvScalarMatrix TEqn
            //     (
            //         fvm::ddt(rho, Cp, T)
            //      == fvm::laplacian(kappa, T)
            //       + fvOptions(rho*Cp, T)
            //     );

            //     fvOptions.constrain(TEqn);

            //     TEqn.solve();

            //     fvOptions.correct(T);
            // }

           {
               fvVectorMatrix DEqn
               (
                   fvm::d2dt2(D)
                ==
                   sig1*fvm::laplacian(2*(mu1*eta + mu2*(1-eta)) + lambda1*eta + lambda2*(1-eta), D, "laplacian(DD,D)")
                  + (sig1/sig2)*divSigmaExp 
                );  
                
                    DEqn -= (sig1)*fvc::div((2*mu1*eta + 2*mu2*(1-eta)) * eta *cEigenStrain + (lambda1*eta + lambda2*(1 -eta)) * I * tr(eta * cEigenStrain) );
            //     //   + rho*fvOptions.d2dt2(D)
                

            //     // if (thermo.thermalStress())
            //     // {
            //     //     DEqn += fvc::grad(threeKalpha*thermo.T());
            //     // }

            //     // fvOptions.constrain(DEqn);

                initialResidual = DEqn.solve().max().initialResidual();

                if (!compactNormalStress)
                {
                    divSigmaExp = fvc::div(DEqn.flux());
                }
            }

            {
                volTensorField gradD = fvc::grad(D);

                // strain = ((gradD - (eta*eta*eta*(6*eta*eta-15*eta+10)) * cEigenStrain)
                //     && symmTensor(1 0 0 0 0 0)) * symmTensor(1 0 0 0 0 0)
                //     + ((gradD - (eta*eta*eta*(6*eta*eta-15*eta+10)) * cEigenStrain)
                //     && symmTensor(0 0 0 1 0 0)) * symmTensor(0 0 0 1 0 0)
                //     + ((gradD - (eta*eta*eta*(6*eta*eta-15*eta+10))) * cEigenStrain)
                //     && symmTensor(0 0 0 0 0 1)) * symmTensor(0 0 0 0 0 1);

                sigmaD = (mu1*eta + mu2*(1-eta))*twoSymm(gradD) + (lambda1*eta + lambda2*(1 - eta))*(I*tr(gradD));

                if (compactNormalStress)
                {
                    divSigmaExp = sig2*fvc::div
                    (
                        sigmaD - (2*(mu1*eta + mu2 *(1-eta)) + (lambda1*eta + lambda2*(1-eta)))*gradD,
                        "div(sigmaD)"
                    );
                }
                else
                {
                    divSigmaExp += sig2*fvc::div(sigmaD);
                }
            }

        } while (initialResidual > convergenceTolerance && ++iCorr < nCorr);

        #include "calculateStress.H"

        //! Writing the results according to keywords in controlDict
    if (runTime.writeTime())
    {
      runTime.write();
    }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
        runTime.write();
     }

    Info<< "End\n" << endl;

    return 0;

}


// ************************************************************************* //
