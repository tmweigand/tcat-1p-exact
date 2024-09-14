#include "TCAT_incompressible.H" // Include the header file

namespace Foam {

MacroscaleIncompressible::MacroscaleIncompressible(

    const volVectorField &_U, const volScalarField &_p,
    const surfaceScalarField &_phi, const int time)
    : U(_U), p(_p), phi(_phi), mesh(U.mesh()),
      normals(IOobject("normals", U.time().timeName(), U.mesh()), mesh,
              dimensionedVector(dimless, Zero)),
      media_patch(
          mesh.boundaryMesh()[mesh.boundaryMesh().findPatchID("media")]) {

  media_label = mesh.boundaryMesh().findPatchID("media");
  normals = -mesh.Sf() / mesh.magSf();

  word file_dir = "tcat/";

  outfile_ptr.reset(
      new OFstream(file_dir + "test_" + std::to_string(time) + ".txt"));
}

void MacroscaleIncompressible::update(const scalar time) {

  singlePhaseTransportModel laminar_transport(U, phi);
  turbulence =
      (incompressible::turbulenceModel::New(U, phi, laminar_transport));

  dimensionedScalar micro_mass_error(
      "micro_mass_error", mass_error().dimensions(),
      returnReduce(mass_error().value(), sumOp<double>()));

  dimensionedVector micro_mom_error(
      "micro_mom_error", momentum_error().dimensions(),
      returnReduce(momentum_error().value(), sumOp<Vector<double>>()));
  momentum_error();

  dimensionedVector time_(
      "time_int", time_int().dimensions(),
      returnReduce(time_int().value(), sumOp<Vector<double>>()));

  dimensionedVector g_p(
      "g_p", grad_p().dimensions(),
      returnReduce(grad_p().value(), sumOp<Vector<double>>()));

  dimensionedVector div_stress_tensor(
      "div_stress_tensor", div_st().dimensions(),
      returnReduce(div_st().value(), sumOp<Vector<double>>()));

  dimensionedVector div_velocity(
      "div_velocity", div_vel().dimensions(),
      returnReduce(div_vel().value(), sumOp<Vector<double>>()));  

  dimensionedVector surface(
      "surface", surface_int().dimensions(),
      returnReduce(surface_int().value(), sumOp<Vector<double>>()));

  dimensionedScalar surface_area(
      "surface_area", surface_area_int().dimensions(),
      returnReduce(surface_area_int().value(), sumOp<double>()));

  dimensionedVector macro_mom_error = time_ + div_velocity + g_p + div_stress_tensor + surface;

  // Output
  outfile_ptr() << "Time: " << time << "\n"
                << "Surface Area: " << surface_area.value() << "\n"
                << "Microscale mass error: " << micro_mass_error.value() << "\n"
                << "Microscale momentum error: " << micro_mom_error.value()
                << "\n"
                << "Macroscale momentum error: " << macro_mom_error.value()
                << "\n"
                << "Macroscale time integral: " << time_.value() << "\n"
                << "Macroscale pressure gradient: " << g_p.value() << "\n"
                << "Macroscale div stress tensor: " << div_stress_tensor.value()
                << "\n"
                << "Macroscale div velocity: " << div_velocity.value()
                << "\n"
                << "Interphase momentum transfer: " << surface.value() << "\n";
}

tmp<volSymmTensorField> MacroscaleIncompressible::get_stress_tensor() {

  tmp<volScalarField> tnuEff = turbulence->nuEff();

  return tmp<volSymmTensorField>::New(tnuEff() * devTwoSymm(fvc::grad(U)));
}

tmp<volVectorField> MacroscaleIncompressible::get_div_stress_tensor() {

  tmp<volVectorField> tU(U);
  tmp<volScalarField> tnuEff = turbulence->nuEff();

  return tmp<volVectorField>::New(
      "divDevRhoReff",
      -fvc::div(tnuEff() * dev2(T(fvc::grad(tU()))),
                "div((nuEff*dev2(T(grad(U)))))") -
          fvc::laplacian(tnuEff(), tU(), "laplacian(nuEff,U)"));
}

dimensionedVector MacroscaleIncompressible::time_int() {
  return (fvc::domainIntegrate(fvc::ddt(U)));
}

dimensionedVector MacroscaleIncompressible::grad_p() {
  return (fvc::domainIntegrate(fvc::grad(p)) -
          surface_integrate(fvc::interpolate(p) * normals));
}

dimensionedVector MacroscaleIncompressible::div_st() {
  return (fvc::domainIntegrate(get_div_stress_tensor()) -
          surface_integrate(fvc::interpolate(get_stress_tensor()) & normals));
}

dimensionedVector MacroscaleIncompressible::div_vel() {
  return (fvc::domainIntegrate(fvc::div(phi,U)) -
          surface_integrate(fvc::interpolate(U*U) & normals));
}



dimensionedVector MacroscaleIncompressible::surface_int() {
  return (surface_integrate(fvc::interpolate(p) * normals) +
          surface_integrate(fvc::interpolate(get_stress_tensor()) & normals));
}

dimensionedVector MacroscaleIncompressible::momentum_error() {
  return (fvc::domainIntegrate(fvc::ddt(U) + fvc::div(phi, U) + fvc::grad(p) +
                               get_div_stress_tensor()));
}

dimensionedScalar MacroscaleIncompressible::mass_error() {
  return (fvc::domainIntegrate(fvc::div(phi)));
}

dimensionedVector
MacroscaleIncompressible::surface_integrate(const surfaceVectorField &data) {

  const surfaceScalarField &magSf = mesh.magSf();
  auto data_dim = data.dimensions();
  data_dim[1] = data_dim[1] + magSf.dimensions()[1];
  vector surface_int(0, 0, 0);
  forAll(media_patch, face) {
    surface_int += data.boundaryField()[media_label][face] *
                   magSf.boundaryField()[media_label][face];
  }

  dimensionedVector dim_surface_int("surface_int", data_dim, surface_int);
  return dim_surface_int;
}

dimensionedScalar MacroscaleIncompressible::surface_area_int() {

  const surfaceScalarField &magSf = mesh.magSf();
  scalar surface_int = 0;
  forAll(media_patch, face) {
    surface_int += magSf.boundaryField()[media_label][face];
  }

  dimensionedScalar dim_surface_int("surface_int", magSf.dimensions(),
                                    surface_int);
  return dim_surface_int;
}

} // namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //