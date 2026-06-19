if (runTime.writeTime())
{
    volSymmTensorField sigma
    (
        IOobject
        (
            "sigma",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rho*sigmaD
    );  
    
    sigma = sigma - (lambda1*eta +lambda2*(1-eta))*rho*I*tr(eta*cEigenStrain)-(rho*(2*(mu1*eta +mu2*(1-eta))*eta*cEigenStrain));
    

    // if (thermo.thermalStress())
    // {
    //     sigma = sigma - I*(threeKalpha*thermo.T());
    // }

    // Calculate Hydrostatic Pressure: -(1/3) * trace(sigma)
    volScalarField hydrostaticPressure
    (
        IOobject
        (
            "hydrostaticPressure",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        -(1.0/3.0)*tr(sigma)
    );

    Info<< "Max hydrostaticPressure = " << max(hydrostaticPressure).value() << endl;



    volScalarField sigmaEq
    (
        IOobject
        (
            "sigmaEq",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt((3.0/2.0)*magSqr(dev(sigma)))
    );

    Info<< "Max sigmaEq = " << max(sigmaEq).value()
        << endl;

    runTime.write();
}
