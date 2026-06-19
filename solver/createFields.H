#include "readTransportProperties.H"
#include "readMechanicalProperties.H"

// IOdictionary thermalProperties
// (
//     IOobject
//     (
//         "thermalProperties",
//         runTime.constant(),
//         mesh,
//         IOobject::MUST_READ_IF_MODIFIED,
//         IOobject::NO_WRITE
//     )
// );

volScalarField eta
(
    IOobject
    (
        "eta",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volScalarField muC
(
    IOobject
    (
        "muC",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
volScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
volScalarField c
(
    IOobject
    (
        "c",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);



// solidDisplacementFoam



Info<< "Reading field D\n" << endl;
volVectorField D
(
    IOobject
    (
        "D",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volTensorField gradD
(
    IOobject
    (
        "gradD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::grad(D)
);



volScalarField etaOld
(
    IOobject
    (
        "etaOld",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

IOdictionary constProperties
    (
        IOobject
        (
            "constProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );


    Info<< "Reading diffusivity DT\n" << endl;

    dimensionedScalar sig
    (
        constProperties.lookup("sig")
    );
  
    dimensionedScalar lambdas
    (
        constProperties.lookup("lambdas")
    );
    dimensionedScalar DTs
    (
        constProperties.lookup("DTs")
    );
    dimensionedSymmTensor cEigenStrain
    (
        constProperties.lookup("cEigenStrain")
    );
    dimensionedScalar mu1
    (
        constProperties.lookup("mu1")
    );
    dimensionedScalar mu2
    (
        constProperties.lookup("mu2")
    );
    dimensionedScalar lambda1
    (
        constProperties.lookup("lambda1")
    );
    dimensionedScalar lambda2
    (
        constProperties.lookup("lambda2")
    );
    dimensionedScalar sig1
    (
        constProperties.lookup("sig1")
    );
    dimensionedScalar sig2
    (
        constProperties.lookup("sig2")
    );



Info<< "Calculating stress field sigmaD\n" << endl;
volSymmTensorField sigmaD
(
    IOobject
    (
        "sigmaD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    (mu1*eta + mu2*(1-eta))*twoSymm(fvc::grad(D)) + (lambda1*eta +lambda2*(1-eta))*(I*tr(fvc::grad(D)))
);

volSymmTensorField deltaSigmaD
(
    IOobject
    (
        "deltaSigmaD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
     ((mu1-mu2)*twoSymm(fvc::grad(D)) + (lambda1-lambda2)*(I*tr(fvc::grad(D))) - 2*(mu1-mu2)*eta*cEigenStrain - (lambda1-lambda2)*(I*tr(eta*cEigenStrain)))
);

volSymmTensorField Sigma
(
    IOobject
    (
        "Sigma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    (2*(mu1*eta +mu2*(1-eta))*(0.5*twoSymm(fvc::grad(D)) - eta*cEigenStrain) + (lambda1*eta +lambda2*(1-eta))*(I*tr(fvc::grad(D) - eta*cEigenStrain)))
);


Info<< "Calculating explicit part of div(sigma) divSigmaExp\n" << endl;
volVectorField divSigmaExp
(
    IOobject
    (
        "divSigmaExp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    sig2*fvc::div(sigmaD)
);

volScalarField effCo
(
    IOobject
    (
        "effCo",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    2.0*(mu1*eta + mu2*(1.0 - eta)) + lambda1*eta + lambda2*(1.0 - eta)
);

if (compactNormalStress)
{
    divSigmaExp -= sig2*fvc::laplacian(effCo, D, "laplacian(DD,D)");
}
else
{
    divSigmaExp -= sig2*fvc::div(effCo*gradD, "div(sigmaD)");
}

mesh.setFluxRequired(D.name());

#include "createFvOptions.H"


