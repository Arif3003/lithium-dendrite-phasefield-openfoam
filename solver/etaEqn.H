           
           
            Random obj(1);
            const scalar randNumber(obj.scalar01());

            Sigma = (2*(mu1*eta+mu2*(1-eta))*(0.5*twoSymm(fvc::grad(D)) - eta*cEigenStrain) + (lambda1*eta +lambda2*(1-eta))*(I*tr(fvc::grad(D) - eta*cEigenStrain)));

            deltaSigmaD = ((mu1-mu2)*twoSymm(fvc::grad(D)) + (lambda1-lambda2)*(I*tr(fvc::grad(D))) - 2*(mu1-mu2)*eta*cEigenStrain - (lambda1-lambda2)*(I*tr(eta*cEigenStrain)));


            volScalarField h_e = 30*eta*eta*(1-eta)*(1-eta);
            volScalarField c_Lip = ((exp((muC-(e_l*R*T))/(R*T)))/(1+exp((muC-(e_l*R*T))/(R*T))))*(1-(eta*eta*eta*(6*eta*eta-15*eta+10)));
            volScalarField c_l= (exp((muC-(e_l*R*T))/(R*T))/(1+exp((muC-(e_l*R*T))/(R*T))));
            volScalarField c_s= (exp((muC-(e_s*R*T))/(R*T))/(1+exp((muC-(e_s*R*T))/(R*T))));
            volScalarField chi= (1/(R*T))*(((exp((muC-(e_l*R*T))/(R*T)))/((1+exp((muC-(e_l*R*T))/(R*T)))*(1+exp((muC-(e_l*R*T))/(R*T)))))*(1-(eta*eta*eta*(6*eta*eta-15*eta+10)))+((C_m_s/C_m_l)*((exp((muC-(e_s*R*T))/(R*T)))/((1+exp((muC-(e_s*R*T))/(R*T)))*(1+exp((muC-(e_s*R*T))/(R*T)))))*(eta*eta*eta*(6*eta*eta-15*eta+10))));
            volScalarField sigma1= S_s*(eta*eta*eta*(6*eta*eta-15*eta+10))+S_l*(1-(eta*eta*eta*(6*eta*eta-15*eta+10)));
            volScalarField Deff = D_s*(eta*eta*eta*(6*eta*eta-15*eta+10))+D_l*(1-(eta*eta*eta*(6*eta*eta-15*eta+10)));
            
            Info<< "Min/max h_e:" << min(h_e).value() << ' '
            << max(h_e).value() << endl;
            Info<< "Min/max c_Lip:" << min(c_Lip).value() << ' '
            << max(c_Lip).value() << endl;
            Info<< "Min/max c_s:" << min(c_s).value() << ' '
            << max(c_s).value() << endl;
            Info<< "Min/max c_l:" << min(c_l).value() << ' '
            << max(c_l).value() << endl;
            Info<< "Min/max chi:" << min(chi).value() << ' '
            << max(chi).value() << endl;
            Info<< "Min/max sigma1:" << min(sigma1).value() << ' '
            << max(sigma1).value() << endl;
                       
           #include "dAdgradEtaMod.H"

            fvScalarMatrix etaEqn
            (
              dimt*fvm::ddt(eta) == -L_s*(-1.0*k*dimx*dimx*fvm::laplacian(ac_01*ac_01,eta) - 1.0*k*dimx*fvc::div(dadgradPhi) + 2*W*eta*(1-eta)*(1-2*eta) - 6*noise_mag*eta*(1.0-eta)*eta*(1-eta)*randNumber)-L_e*h_e*(exp((1-alpha)*n*(F*phi)/(R*T))-((C_m_l*c_Lip)/C_0)*exp(-(alpha)*n*(F*phi)/(R*T))) 
	     );
            etaEqn += (L_s)*(1)*(0.5*(deltaSigmaD && (symm(fvc::grad(D))- eta*cEigenStrain))-(Sigma && cEigenStrain));
            etaEqn.solve();
            volVectorField n1=dimx*fvc::grad(eta)/(1E-20+mag(dimx*fvc::grad(eta)));
            
            fvScalarMatrix muCEqn
             (
               (chi)*dimt*fvm::ddt(muC) == dimx*dimx*fvm::laplacian(((D1*c_Lip)/(R*T)),muC)+dimx*dimx*fvc::laplacian(((D1*c_Lip*n*F)/(R*T)),phi) - dimt*(h_e)*(c_s*(C_m_s/C_m_l)-c_l)*fvc::ddt(eta)             
	     );            
             muCEqn.solve();
             Info<< "Min/max eta:" << min(eta).value() << ' '
            << max(eta).value() << endl;
             Info<< "Min/max muC:" << min(muC).value() << ' '
            << max(muC).value() << endl;
                       
            solve
            (
               (1/(dx*dx*1e-12))*fvm::laplacian(sigma1,phi) == (n*F*C_m_s_r/(4000))*fvc::ddt(eta)                
            );
            Info<< "Min/max phi:" << min(phi).value() << ' '
            << max(phi).value() << endl;
// ************************************************************************* //
