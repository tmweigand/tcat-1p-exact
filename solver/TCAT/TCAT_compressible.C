#include "TCAT_compressible.H" // Include the header file

namespace Foam {

MacroscaleCompressible::MacroscaleCompressible(
    word file_out, const volVectorField &_U, const volScalarField &_p,
    const volScalarField &_p_rgh, const volScalarField &_rho,
    const surfaceScalarField &_phi, const volScalarField &_chem_potential, 
    const volScalarField &_grav_potential, const meshObjects::gravity& _g, const volScalarField &_gh,
    const surfaceScalarField &_ghf, dimensionedScalar _v_in,
    dimensionedScalar _rho0, dimensionedScalar _beta, dimensionedScalar _domain_volume,
    dimensionedScalar _mu, const int time)
    : U(_U), p(_p), p_rgh(_p_rgh), rho(_rho), phi(_phi), chem_potential(_chem_potential), grav_potential(_grav_potential), g(_g), gh(_gh), ghf(_ghf),
      v_in(_v_in), rho0(_rho0), beta(_beta), domain_volume(_domain_volume), mu(_mu),
      mesh(U.mesh()),
      normals(IOobject("normals", U.time().timeName(), U.mesh()), mesh,
              dimensionedVector(dimless, Zero)),
      media_patch(
          mesh.boundaryMesh()[mesh.boundaryMesh().findPatchID("media")]) {

  media_label = mesh.boundaryMesh().findPatchID("media");
  normals = mesh.Sf() / mesh.magSf(); // n_w - outward normal of w_phase

  word file_dir = "tcat/tcat_out_";

  outfile_ptr.reset(
      new OFstream(file_dir + file_out + "_" + std::to_string(time) + ".txt"));
}

void MacroscaleCompressible::update(const scalar time) {

  singlePhaseTransportModel laminar_transport(U, phi);

  turbulence =
      (incompressible::turbulenceModel::New(U, phi, laminar_transport));

  dimensionedScalar w_vol("w_vol", dimVolume, gSum(mesh.V()));

  dimensionedScalar ws_area(get_surface_area());


  // Macroscale Variables 
  dimensionedScalar m_density("m_density", macro_density() / w_vol);

  dimensionedVector m_velocity("m_velocity", macro_velocity());

  dimensionedScalar m_pressure("m_pressure", macro_pressure() / w_vol);

  dimensionedVector m_gravity("m_gravity", macro_gravity());

  dimensionedScalar m_chem_potential ("m_chem_potential", macro_chem_potential());

  dimensionedScalar m_grav_potential ("m_grav_potential", macro_grav_potential());

  // Momentum  Terms 
  dimensionedVector time_integral(time_int() / domain_volume);

  dimensionedVector div_velocity(div_vel() / domain_volume);

  dimensionedVector grad_pressure( grad_p()/domain_volume + e_w*m_density*m_gravity );

  dimensionedVector div_stress_tensor(-grad_pressure + div_st()/domain_volume);

  dimensionedVector surface(surface_int() / domain_volume);

  // dimensionedVector grad_chem("grad_chem", grad_chem_potential() / w_vol );

  // dimensionedVector grad_potential(grad_psi() /w_vol);

  // dimensionedVector grad_ew_rho("grad_e_rho", grad_e_rho() );  

  dimensionedVector e_rho_grad_chem("e_rho_grad_chem", grad_e_rho_chem() - m_chem_potential * grad_e_rho() );  

  dimensionedVector grad_chem("grad_chem", grad_chem_potential() );

  dimensionedVector grad_potential(grad_psi() );

  dimensionedTensor m_stress_tensor(
      "m_stress_tensor", macro_stress_tensor(m_pressure, m_velocity) / w_vol);

  dimensionedTensor m_grad_u("grad_u", grad_u() / w_vol);

  dimensionedTensor m_strain_tensor("grad_u", macro_strain_tensor(m_grad_u));

  // Errors 
  dimensionedScalar micro_mass_error(mass_error()/domain_volume);

  dimensionedVector micro_mom_error(momentum_error()/domain_volume);

  dimensionedVector macro_mom_error = (time_integral + div_velocity - div_stress_tensor - e_w*m_density*m_gravity - surface);

  e_w = w_vol / domain_volume;
  dimensionedScalar e_ws("e_ws", dimless,
                         ws_area.value() / domain_volume.value());
  dimensionedScalar sauder_mean("sauder_mean", dimLength,
                                6 * (1 - e_w.value()) / e_ws.value());

  dimensionedScalar Re("Re", reynolds(sauder_mean, m_density, mag(m_velocity)));

  tensor Ei(1, 0, 0, 0, 1, 0, 0, 0, 1);
  dimensionedScalar approx_1 =
      (e_w * m_stress_tensor + e_w * m_pressure * Ei) && m_strain_tensor;
  dimensionedScalar fit = -(e_rho_grad_chem[0] + e_w * m_density*grad_potential[0]) /
                          m_velocity[0];
  dimensionedScalar fit_2 = -(e_w*grad_pressure[0] + e_w*m_density*m_gravity[0])/m_velocity[0];

  dimensionedScalar K = e_w * e_w * m_density / fit;
  dimensionedScalar K2 = e_w * e_w * m_density / fit_2;

  // Output
  outfile_ptr() << "Time," << time << "\n"
                << "Reynolds Number," << Re.value() << "," << Re.dimensions() << '\n'
                << "Fit," << fit.value() << "," << fit.dimensions() << '\n'
                << "K," << K.value() << "," << K.dimensions() << '\n'
                << "TCAT Approx 1," << approx_1.value() << ","
                << approx_1.dimensions() << "\n"
                << "e_w," << e_w.value() << "," << e_w.dimensions() << "\n"
                << "e_ws," << e_ws.value() << "," << e_ws.dimensions() << "\n"
                << "sauder mean diameter," << sauder_mean.value() << ","
                << sauder_mean.dimensions() << "\n"
                << "Domain Volume," << domain_volume.value() << ","
                << domain_volume.dimensions() << "\n"
                << "w_vol," << w_vol.value() << "," << w_vol.dimensions() << '\n'
                << "Macroscale rho," << m_density.value() << ","
                << m_density.dimensions() << "\n"
                << "Macroscale U," << m_velocity.value() << ","
                << m_velocity.dimensions() << "\n"
                << "Macroscale p," << m_pressure.value() << ","
                << m_pressure.dimensions() << "\n"
                << "Macroscale t," << m_stress_tensor.value() << ","
                << m_stress_tensor.dimensions() << "\n"
                << "Macroscale d," << m_strain_tensor.value() << ","
                << m_strain_tensor.dimensions() << "\n"
                << "Macroscale mu," << m_chem_potential.value() << ","
                << m_chem_potential.dimensions() << "\n"
                << "Macroscale psi," << m_grav_potential.value() << ","
                << m_grav_potential.dimensions() << "\n"
                << "Microscale mass error," << micro_mass_error.value() << ","
                << micro_mass_error.dimensions() << "\n"
                << "Microscale momentum error," << micro_mom_error.value()
                << "," << micro_mom_error.dimensions() << "\n"
                << "Macroscale momentum error," << macro_mom_error.value()
                << "," << macro_mom_error.dimensions() << "\n"
                << "Macroscale ddt(rho*U)," << time_integral.value() << ","
                << time_integral.dimensions() << "\n"
                << "Macroscale div(rho*U*U)," << div_velocity.value() << ","
                << div_velocity.dimensions() << "\n"
                << "Macroscale rho*g," << (m_density * m_gravity).value()
                << "," << (m_density * m_gravity).dimensions() << "\n"
                << "Macroscale div(t)," << div_stress_tensor.value() << ","
                << div_stress_tensor.dimensions() << "\n"
                << "Macroscale T," << surface.value() << ","
                << surface.dimensions() << "\n"
                << "Macroscale grad(pressure),"
                << grad_pressure.value() << ","
                << grad_pressure.dimensions() << "\n"
                << "Macroscale e_w*rho*grad(chem potential),"
                << (e_w*m_density*grad_chem).value() << ","
                << (e_w*m_density*grad_chem).dimensions() << "\n"
                << "Macroscale ew*rho*grad(potential)," << (e_w*m_density*grad_potential).value()
                << "," << (e_w*m_density*grad_potential).dimensions() << "\n"
                << "Macroscale Potential Difference," << (e_w*m_density*grad_potential + m_density*grad_potential).value()
                << "," << (m_density*grad_potential + m_density*grad_potential).dimensions() << "\n"
                << "ew*pw*del(mu): " << e_rho_grad_chem.value() << "," << e_rho_grad_chem.dimensions() << "\n"
                ;

  
}

tmp<volSymmTensorField> MacroscaleCompressible::get_stress_tensor() 
{
  return tmp<volSymmTensorField>::New(mu * devTwoSymm(fvc::grad(U)));
}

tmp<volVectorField> MacroscaleCompressible::get_div_stress_tensor() 
{
  tmp<volVectorField> tU(U);
  return tmp<volVectorField>::New(
      "divDevRhoReff",
      fvc::div(mu * dev2(T(fvc::grad(tU()))), "div((mu*dev2(T(grad(U)))))") +
          fvc::laplacian(mu, tU(), "laplacian(mu,U)"));
}

dimensionedVector MacroscaleCompressible::macro_gravity() {
  return (fvc::domainIntegrate(rho * g)/fvc::domainIntegrate(rho));
}

dimensionedVector MacroscaleCompressible::macro_velocity() {
  return (fvc::domainIntegrate(rho * U) / fvc::domainIntegrate(rho));
}

dimensionedScalar MacroscaleCompressible::macro_pressure() {
  return (fvc::domainIntegrate(p));
}

dimensionedScalar MacroscaleCompressible::macro_density() {
  return (fvc::domainIntegrate(rho));
}

dimensionedScalar MacroscaleCompressible::macro_chem_potential(){
  return (fvc::domainIntegrate(rho*chem_potential)/fvc::domainIntegrate(rho));
}

dimensionedScalar MacroscaleCompressible::macro_grav_potential(){
  return (fvc::domainIntegrate(rho*grav_potential)/fvc::domainIntegrate(rho));
}

dimensionedTensor
MacroscaleCompressible::macro_stress_tensor(dimensionedScalar macro_pressure,
                                            dimensionedVector macro_velocity) {
  tensor Ei(1, 0, 0, 0, 1, 0, 0, 0, 1);

  return (
      fvc::domainIntegrate(-p * Ei + get_stress_tensor() -
                           rho * (U - macro_velocity) * (U - macro_velocity)));
}

dimensionedTensor
MacroscaleCompressible::macro_strain_tensor(dimensionedTensor macro_grad_u) {

  return (0.5 * (macro_grad_u + macro_grad_u.T()));
}

dimensionedScalar
MacroscaleCompressible::reynolds(dimensionedScalar sauder_mean,
                                 dimensionedScalar rho, dimensionedScalar u) {
  return ((rho * sauder_mean * u) / mu);
}


dimensionedVector MacroscaleCompressible::time_int() {
  return (fvc::domainIntegrate(fvc::ddt(rho, U)));
}

dimensionedVector MacroscaleCompressible::grad_e_rho(){
  return ((fvc::domainIntegrate(fvc::grad(rho)) -
          surface_integrate(fvc::interpolate(rho) * normals))/domain_volume);
}

dimensionedVector MacroscaleCompressible::grad_e_rho_chem(){
  return ((fvc::domainIntegrate(fvc::grad(rho*chem_potential)) -
          surface_integrate(fvc::interpolate(rho*chem_potential) * normals))/domain_volume);
}

dimensionedTensor MacroscaleCompressible::grad_u() {
  return (fvc::domainIntegrate(fvc::grad(U)) -
          surface_integrate(fvc::interpolate(U) * normals));
}

dimensionedVector MacroscaleCompressible::grad_p() {
  // Reconstruct  form p_rgh rather than p to mimic what OF solves
  // return (fvc::domainIntegrate(fvc::grad(p)) -
  // surface_integrate(fvc::interpolate(p) * normals));
  tensor Ei(1, 0, 0, 0, 1, 0, 0, 0, 1);
  return (fvc::domainIntegrate(fvc::reconstruct(
              (ghf * fvc::snGrad(rho) + fvc::snGrad(p_rgh)) * mesh.magSf())) -
          surface_integrate(fvc::interpolate(p*Ei) & normals));
}

dimensionedVector MacroscaleCompressible::grad_p_new() {
  // Reconstruct  form p_rgh rather than p to mimic what OF solves
  // return (fvc::domainIntegrate(fvc::grad(p)) -
  // surface_integrate(fvc::interpolate(p) * normals));
  tensor Ei(1, 0, 0, 0, 1, 0, 0, 0, 1);

  dimensionedVector scale("scale", dimensionSet(0, -1, 0, 0, 0, 0, 0),
                          vector(1, 0, 0));
  return (fvc::domainIntegrate(rho*gh*scale) + fvc::domainIntegrate(fvc::reconstruct(
              (ghf * fvc::snGrad(rho) + fvc::snGrad(p_rgh)) * mesh.magSf())) -
          surface_integrate(fvc::interpolate(p*Ei) & normals));
}


dimensionedVector MacroscaleCompressible::div_st() {
  tensor Ei(1, 0, 0, 0, 1, 0, 0, 0, 1);

  return (fvc::domainIntegrate(get_div_stress_tensor()) -
          surface_integrate(fvc::interpolate(get_stress_tensor()) & normals));
}

dimensionedVector MacroscaleCompressible::div_vel() {
  return (fvc::domainIntegrate(fvc::div(phi, U)) -
          surface_integrate(fvc::interpolate(rho * U * U) & normals));
}

dimensionedVector MacroscaleCompressible::surface_int() {
  tensor Ei(1, 0, 0, 0, 1, 0, 0, 0, 1);
  return (-surface_integrate(fvc::interpolate(p*Ei) & normals) 
          + surface_integrate(fvc::interpolate(get_stress_tensor()) & normals));
}


dimensionedVector MacroscaleCompressible::grad_chem_potential() {

  // return (fvc::domainIntegrate(fvc::grad(chem_potential)) -
  //         surface_integrate(fvc::interpolate(chem_potential) * normals));
    return ( (fvc::domainIntegrate(fvc::grad(rho*chem_potential)) -
          surface_integrate(fvc::interpolate(rho*chem_potential) * normals)) / fvc::domainIntegrate(rho) );
}


dimensionedVector MacroscaleCompressible::grad_psi() {
  // return (fvc::domainIntegrate(fvc::grad(grav_potential)) -
  //         surface_integrate(fvc::interpolate(grav_potential) * normals));
    return ( (fvc::domainIntegrate(fvc::grad(rho*grav_potential)) -
          surface_integrate(fvc::interpolate(rho*grav_potential) * normals)) / fvc::domainIntegrate(rho) );
}

// dimensionedVector MacroscaleCompressible::grad_psi() {

//   return (fvc::domainIntegrate(fvc::grad(rho * gh)) -
//           surface_integrate(fvc::interpolate(rho) * ghf * normals));
// }

dimensionedVector MacroscaleCompressible::grad_psi_2() {

  dimensionedVector scale("scale", dimensionSet(0, -1, 0, 0, 0, 0, 0),
                          vector(1, 0, 0));

  return (fvc::domainIntegrate(rho * gh * scale) +
          fvc::domainIntegrate(
              fvc::reconstruct((ghf * fvc::snGrad(rho)) * mesh.magSf())) -
          surface_integrate(fvc::interpolate(rho) * ghf * normals));
}

dimensionedVector MacroscaleCompressible::momentum_error() {
  return (fvc::domainIntegrate(
      fvc::ddt(rho, U) + fvc::div(phi, U) - get_div_stress_tensor() +
      fvc::reconstruct((ghf * fvc::snGrad(rho) + fvc::snGrad(p_rgh)) *
                       mesh.magSf())));
}

dimensionedScalar MacroscaleCompressible::mass_error() {
  return (fvc::domainIntegrate(fvc::ddt(rho) + fvc::div(phi)));
}



dimensionedVector
MacroscaleCompressible::surface_integrate(const surfaceVectorField &data) {

  const surfaceScalarField &magSf = mesh.magSf();
  auto data_dim = data.dimensions();
  data_dim[1] = data_dim[1] + magSf.dimensions()[1];
  vector surface_int(0, 0, 0);
  forAll(media_patch, face) {
    surface_int += data.boundaryField()[media_label][face] *
                   magSf.boundaryField()[media_label][face];
  }

  dimensionedVector dim_surface_int(
      "surface_int", data_dim,
      returnReduce(surface_int, sumOp<Vector<double>>()));

  return dim_surface_int;
}

dimensionedTensor
MacroscaleCompressible::surface_integrate(const surfaceTensorField &data) {

  const surfaceScalarField &magSf = mesh.magSf();
  auto data_dim = data.dimensions();
  data_dim[1] = data_dim[1] + magSf.dimensions()[1];
  tensor surface_int(0, 0, 0, 0, 0, 0, 0, 0, 0);
  forAll(media_patch, face) {
    surface_int += data.boundaryField()[media_label][face] *
                   magSf.boundaryField()[media_label][face];
  }

  dimensionedTensor dim_surface_int(
      "surface_int", data_dim,
      returnReduce(surface_int, sumOp<Tensor<double>>()));
  return dim_surface_int;
}

dimensionedScalar MacroscaleCompressible::get_surface_area() {

  const surfaceScalarField &magSf = mesh.magSf();
  scalar surface_int = 0;
  forAll(media_patch, face) {
    surface_int += magSf.boundaryField()[media_label][face];
  }

  dimensionedScalar dim_surface_int("surface_int", magSf.dimensions(),
                                    returnReduce(surface_int, sumOp<double>()));
  return dim_surface_int;
}

} // namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //