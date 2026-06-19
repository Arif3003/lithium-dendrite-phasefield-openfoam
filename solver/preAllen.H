fvScalarMatrix etaEqn
            (
                dimt*fvm::ddt(eta) == -L_s*(-dimx*dimx*(1)*fvm::laplacian(k,eta) + 2*W*eta*(1-eta)*(1-2*eta))



            //   dimt*fvm::ddt(eta) == -L_s*(-2.0*k*dimx*dimx*fvm::laplacian(ac_01*ac_01,eta) + 2.0*k*dimx*fvc::div(dadgradphi) + 2*W*eta*(1-eta)*(1-2*eta) - 6*noise_mag*phi*(1.0-phi)*phi*(1-phi)*randNumber)-L_e*h_e*(exp((1-alpha)*n*(F*phi)/(R*T))-((C_m_l*c_Lip)/C_0)*exp(-(alpha)*n*(F*phi)/(R*T))) 
	     );


etaEqn.solve();

Info<< "Min/max eta (Smoothing): " << min(eta).value() << " " << max(eta).value() << endl;