

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

word file_out(transportProperties.lookup("file_out"));

dimensionedScalar mu(transportProperties.lookup("mu"));
dimensionedScalar rho0(transportProperties.lookup("rho0"));
dimensionedScalar beta(transportProperties.lookup("beta"));
dimensionedScalar p0(transportProperties.lookup("p0"));

dimensionedScalar domain_volume(transportProperties.lookup("domain_volume"));

dimensionedScalar v_in(transportProperties.lookup("v_in"));


Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading rho\n" << endl;
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rho0*exp(beta*(p-p0))
);


volScalarField psi
(
    IOobject
    (
        "psi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho0*beta*exp(beta*(p-p0))
);

#include "compressibleCreatePhi.H"
pressureControl pressureControl(p, rho, pimple.dict(), false);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Force p_rgh to be consistent with p
p_rgh = p - rho*gh;

mesh.setFluxRequired(p_rgh.name());

label pRefCell = 0;
scalar pRefValue = 0.0;

if (p_rgh.needReference())
{
    setRefCell
    (
        p,
        p_rgh,
        pimple.dict(),
        pRefCell,
        pRefValue
    );

    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
}

dimensionedScalar initialMass("initialMass", fvc::domainIntegrate(rho));

IOobject dpdtHeader
(
    "dpdt",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
);

volScalarField dpdt(dpdtHeader, fvc::ddt(p));

#include "createMRF.H"
#include "createFvOptions.H"


const dimensionedScalar rhoMax("rhoMax", dimDensity, GREAT, pimple.dict());
const dimensionedScalar rhoMin("rhoMin", dimDensity, Zero, pimple.dict());



volScalarField chem_potential
(
    IOobject
    (
        "chem_potential",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    1/(rho0*beta)*(exp(-beta*p) - exp(-beta*p0))
);

volScalarField grav_potential
(
    IOobject
    (
        "grav_potential",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    gh
);