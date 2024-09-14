tensor Ei(1,0,0,0,1,0,0,0,1);

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


volTensorField stress_t
(
    IOobject
    (
        "stress_t",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p*Ei //Just to initiate
);


volTensorField gradU
(
    IOobject
    (
        "gradU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::grad(U)
);

volVectorField gradp
(
    IOobject
    (
        "gradp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::grad(p)
);


volTensorField strain_tensor
(
    IOobject
    (
        "strain_tensor",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    gradU
);

 surfaceVectorField ten_dot_norm
(
    IOobject
    (
        "ten_dot_norm",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector(dimensionSet(1, -1, -2, 0, 0, 0, 0), Zero)
);

surfaceVectorField normals
(
   IOobject
   (
       "normals",
       runTime.timeName(),
       mesh,
       IOobject::NO_READ,
       IOobject::AUTO_WRITE
   ),
   mesh,
   dimensionedVector(dimensionSet(0, 0, 0, 0, 0, 0, 0), Zero)
);

surfaceVectorField stress_term_dev
(
    IOobject
    (
        "stress_term_dev",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector(dimensionSet(1, -1, -2, 0, 0, 0, 0), Zero)
);


surfaceVectorField stress_term_dev_macro
(
    IOobject
    (
        "stress_term_dev_macro",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector(dimensionSet(1, -1, -2, 0, 0, 0, 0), Zero)
);

volScalarField div_u
(
    IOobject
    (
        "div_u",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::div(U) //Just to initiate
);

volVectorField div_stress_term_dev
(
    IOobject
    (
        "div_stress_term_dev",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::div(stress_t) //Just to initiate
);

volVectorField div_stress_term_dev_macro
(
    IOobject
    (
        "div_stress_term_dev_macro",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::div(stress_t) //Just to initiate
);

volVectorField ddt_u
(
    IOobject
    (
        "ddt_u",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::ddt(U) //Just to initiate
);


surfaceVectorField mass_flux_surface
(
    IOobject
    (
        "mass_flux_surface",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector(dimensionSet(1, -2, -1, 0, 0, 0, 0), Zero)
);


volVectorField div_stress_t
(
    IOobject
    (
        "div_stress_t",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::div(stress_t) //Just to initiate
);

volVectorField grad_p
(
    IOobject
    (
        "grad_p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::grad(p) //Just to initiate
);

surfaceVectorField surface_p
(
    IOobject
    (
        "surface_p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::interpolate(p)*normals //Just to initiate
);

volVectorField adv
(
    IOobject
    (
        "adv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::div(phi, U) //Just to initiate
);

surfaceVectorField surface_stress_t
(
    IOobject
    (
        "surface_stress_t",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::interpolate(stress_t) & normals //Just to initiate
);

volVectorField momentErr
(
    IOobject
    (
        "momentErr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::grad(p) //Just to initiate
);
