#include "TCAT_compressible.H" // Include the header file

namespace Foam
{

MacroscaleCompressible::MacroscaleCompressible(bool local,
                                               word file_out,
                                               const fvMesh& mesh,
                                               const int time,
                                               dimensionedScalar _mu)
    : local(local),
      file_out(file_out),
      mesh(mesh),
      media_label(mesh.boundaryMesh().findPatchID("media")),
      media_patch(
          mesh.boundaryMesh()[mesh.boundaryMesh().findPatchID("media")]),
      mu(_mu),
      n_w(mesh.Sf() / mesh.magSf()),
      averaging_region(mesh.points(), !local),
      domain_volume("domain_volume", dimVolume, averaging_region.volume()),
      domain_center(averaging_region.centre()),
      w_volume( get_volume() ),
      ws_area(get_surface_area()),
      e_w(w_volume / domain_volume),
      e_ws(ws_area / domain_volume),
      sauder_mean(6 * (1 - e_w) / e_w)
{
  
    if (local){
        word file_dir = "tcat/local/";
        // mkDir(file_dir);
    }

}

void
MacroscaleCompressible::update(const int iter,
                               scalar time,
                               const volVectorField& U,
                               const volScalarField& p,
                               const volScalarField& p_rgh,
                               const volScalarField& rho,
                               const surfaceScalarField& phi,
                               const volScalarField& chem_potential,
                               const volScalarField& grav_potential,
                               const meshObjects::gravity& g,
                               const volScalarField& gh,
                               const surfaceScalarField& ghf)
{
    int proc_id = Pstream::myProcNo();

    // Output - needs to be improved
    if (local)
    {
        word file_dir = "tcat/local/";
        // word file_dir = "tcat/local/"+std::to_string(iter)+"/";
        // mkDir(file_dir);
        outfile_ptr.reset(new OFstream(file_dir + file_out + "_" +
                                       std::to_string(iter) + "_" +
                                       std::to_string(proc_id) + ".txt"));

        vtk_ptr.reset(new OFstream(file_dir + file_out + "_VTK_" +
                                   std::to_string(iter) + "_" +
                                   std::to_string(proc_id) + ".vtk"));

        vtk_write_domain();
    }
    else
    {
        word file_dir = "tcat/";
        outfile_ptr.reset(new OFstream(file_dir + file_out + "_" +
                                       std::to_string(iter) + "_global.txt"));
    }
    // Helpful Variables
    tensor Ei(1, 0, 0, 0, 1, 0, 0, 0, 1);
    volTensorField stress_tensor(-p * Ei + get_stress_tensor(mu, U));
    tmp<volVectorField> micro_grad_p(
        rho * g +
        fvc::reconstruct((ghf * fvc::snGrad(rho) + fvc::snGrad(p_rgh)) *
                         mesh.magSf()));

    // Macroscale Variables
    dimensionedScalar m_density = average(rho, w_volume);
    dimensionedScalar m_pressure = average(p, w_volume);

    // Macroscale Variables - Density Weighted
    dimensionedVector m_velocity = average(U, rho);
    dimensionedVector m_gravity = average(g, rho);
    dimensionedScalar m_chem_potential = average(chem_potential, rho);
    dimensionedScalar m_grav_potential = average(grav_potential, rho);

    if (local)
    {
        vtk_write_data("density", m_density);
        vtk_write_data("velocity", m_velocity);
    }

    // Macroscale Stress Tensor
    tmp<volTensorField> stress_tensor_dev =
        stress_tensor - rho * (U - m_velocity) * (U - m_velocity);
    dimensionedTensor m_stress_tensor =
        average(stress_tensor_dev.ref(), w_volume);

    // Mass Terms //
    dimensionedScalar ddt_rho = time_integral(rho);
    dimensionedScalar div_e_rho_v =
        divergence(fvc::div(phi), fvc::interpolate(rho * U));
    dimensionedScalar mass_error = ddt_rho + div_e_rho_v;

    ////////////////////
    // Momentum Terms //
    ////////////////////
    dimensionedVector ddt_rhoU = time_integral(rho, U);
    dimensionedVector div_e_rho_v_v =
        divergence(fvc::div(phi, U), fvc::interpolate(rho * U * U));
    dimensionedVector grad_pressure =
        gradient(micro_grad_p, fvc::interpolate(p));
    dimensionedVector bad_grad_pressure =
        gradient(fvc::grad(p), fvc::interpolate(p));

    dimensionedVector div_stress_tensor =
        -grad_pressure + divergence(get_div_stress_tensor(mu, U),
                                    fvc::interpolate(get_stress_tensor(mu, U)));
    dimensionedVector surface =
        surface_integrate((fvc::interpolate(stress_tensor) & n_w).cref()) /
        domain_volume;

    dimensionedVector mom_error =
        (ddt_rhoU + div_e_rho_v_v - div_stress_tensor -
         e_w * m_density * m_gravity - surface);
    ////////////////////

    ////////////////////
    // Gradient Terms //
    ////////////////////
    dimensionedVector grad_ew_rho =
        gradient(fvc::grad(rho), fvc::interpolate(rho));

    // Gradient Terms - Density Weight so product rule //
    dimensionedVector grad_e_rho_chem =
        gradient( (micro_grad_p + chem_potential*fvc::grad(rho)),
                 fvc::interpolate(rho * chem_potential));
    dimensionedVector e_rho_grad_chem =
        grad_e_rho_chem - (m_chem_potential * grad_ew_rho);


    dimensionedVector grad_e_rho_grav =
        gradient(-rho*g,
                 fvc::interpolate(rho * grav_potential));
    dimensionedVector e_rho_grad_grav =
        grad_e_rho_grav - (m_grav_potential * grad_ew_rho);

    

    dimensionedTensor grad_e_u =
        gradient(fvc::grad(rho * U), fvc::interpolate(rho * U));
    dimensionedTensor e_rho_grad_U = grad_e_u - (m_velocity * grad_ew_rho);
    dimensionedTensor grad_u = e_rho_grad_U / (e_w * m_density);
    dimensionedTensor m_strain_tensor = 0.5 * (grad_u + grad_u.T());
    ////////////////////

    // TCAT Approximations
    dimensionedScalar Re = reynolds(sauder_mean, m_density, mag(m_velocity));

    dimensionedScalar approx_1 =
        (e_w * m_stress_tensor + e_w * m_pressure * Ei) && m_strain_tensor;
    dimensionedScalar fit =
        -(e_rho_grad_chem[0] + e_rho_grad_grav[0]) / m_velocity[0];

    dimensionedScalar K = e_w * e_w * m_density / fit;

    // // Output
    outfile_ptr() << "Time," << time << "\n"
                  << "domain_center," << domain_center << "\n"
                  << "Reynolds Number," << Re.value() << "," << Re.dimensions()
                  << '\n'
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
                  << "w_volume," << w_volume.value() << ","
                  << w_volume.dimensions() << '\n'
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
                  << "Mass error," << mass_error.value() << ","
                  << mass_error.dimensions() << "\n"
                  << "Momentum error," << mom_error.value() << ","
                  << mom_error.dimensions() << "\n"
                  << "Macroscale ddt(rho*U)," << ddt_rhoU.value() << ","
                  << ddt_rhoU.dimensions() << "\n"
                  << "Macroscale div(rho*U*U)," << div_e_rho_v_v.value() << ","
                  << div_e_rho_v_v.dimensions() << "\n"
                  << "Macroscale rho*g," << (m_density * m_gravity).value()
                  << "," << (m_density * m_gravity).dimensions() << "\n"
                  << "Macroscale div(t)," << div_stress_tensor.value() << ","
                  << div_stress_tensor.dimensions() << "\n"
                  << "Macroscale T," << surface.value() << ","
                  << surface.dimensions() << "\n"
                  << "Macroscale grad(pressure)," << grad_pressure.value()
                  << "," << grad_pressure.dimensions() << "\n"
                  << "BAD Macroscale grad(pressure),"
                  << bad_grad_pressure.value() << ","
                  << bad_grad_pressure.dimensions() << "\n"
                  << "Macroscale e_w*rho*grad(chem potential),"
                  << (e_rho_grad_chem).value() << ","
                  << (e_rho_grad_chem).dimensions() << "\n"
                  << "Macroscale e_w*rho*grad(potential),"
                  << (e_rho_grad_grav).value() << ","
                  << (e_rho_grad_grav).dimensions() << "\n";
}

tmp<volSymmTensorField>
MacroscaleCompressible::get_stress_tensor(dimensionedScalar mu,
                                          const volVectorField& U)
{
    return tmp<volSymmTensorField>::New(mu * devTwoSymm(fvc::grad(U)));
}

tmp<volVectorField>
MacroscaleCompressible::get_div_stress_tensor(dimensionedScalar mu,
                                              const volVectorField& U)
{
    tmp<volVectorField> tU(U);
    return tmp<volVectorField>::New(
        "divDevRhoReff",
        fvc::div(mu * dev2(T(fvc::grad(tU()))), "div((mu*dev2(T(grad(U)))))") +
            fvc::laplacian(mu, tU(), "laplacian(mu,U)"));
}

dimensionedScalar
MacroscaleCompressible::reynolds(dimensionedScalar sauder_mean,
                                 dimensionedScalar rho,
                                 dimensionedScalar u)
{
    return ((rho * sauder_mean * u) / mu);
}

dimensionedScalar
MacroscaleCompressible::time_integral(const volScalarField& arg1)
{
    return (volume_integrate(fvc::ddt(arg1).cref()) / domain_volume);
}

dimensionedVector
MacroscaleCompressible::time_integral(const volScalarField& arg1,
                                      const volVectorField& arg2)
{
    return (volume_integrate(fvc::ddt(arg1, arg2).cref()) / domain_volume);
}

template <typename T1, typename T2>
dimensioned<T1>
MacroscaleCompressible::divergence(
    tmp<GeometricField<T1, fvPatchField, volMesh> > div_term,
    tmp<GeometricField<T2, fvsPatchField, surfaceMesh> > surface_term)
{
    return (fvc::domainIntegrate(div_term) -
            surface_integrate((surface_term & n_w).cref())) /
           domain_volume;
}

template <typename T>
dimensioned<T>
MacroscaleCompressible::average(GeometricField<T, fvPatchField, volMesh> data,
                                dimensionedScalar volume)
{
    return (volume_integrate(data) / volume);
}

template <typename T>
dimensioned<T>
MacroscaleCompressible::average(GeometricField<T, fvPatchField, volMesh> data,
                                volScalarField weight)
{
    return (volume_integrate((data * weight).cref()) /
            volume_integrate(weight));
}

dimensionedVector
MacroscaleCompressible::average(const meshObjects::gravity data,
                                volScalarField weight)
{
    return (volume_integrate((data * weight).cref()) /
            volume_integrate(weight));
}

template <typename T>
dimensioned<T>
MacroscaleCompressible::gradient(GeometricField<T, fvPatchField, volMesh> data)
{
    dimensioned<T> _grad(dimVolume * data.dimensions, Zero);
    if (local)
    {
        dimensioned<T> vf_sum(_grad.dimensions(),
                              sum(fvc::volumeIntegrate(data)));
        _grad = vf_sum - surface_integrate((data * n_w).cref());
    }
    else
    {
        _grad =
            fvc::domainIntegrate(data) - surface_integrate((data * n_w).cref());
    }
    return (_grad / domain_volume);
}

template <typename T1, typename T2>
dimensioned<T1>
MacroscaleCompressible::gradient(
    tmp<GeometricField<T1, fvPatchField, volMesh> > grad,
    tmp<GeometricField<T2, fvsPatchField, surfaceMesh> > surface_term)
{
    dimensioned<T1> _grad(dimVolume * grad.ref().dimensions(), Zero);

    if (local)
    {
        dimensioned<T1> vf_sum(_grad.dimensions(),
                               sum(fvc::volumeIntegrate(grad)));
        _grad = vf_sum - surface_integrate((surface_term * n_w).cref());
    }
    else
    {
        _grad = fvc::domainIntegrate(grad) -
                surface_integrate((surface_term * n_w).cref());
    }
    return (_grad / domain_volume);
}

template <typename T>
dimensioned<T>
MacroscaleCompressible::surface_integrate(
    const GeometricField<T, fvsPatchField, surfaceMesh> data)
{
    const surfaceScalarField& magSf = mesh.magSf();
    auto data_dim = data.dimensions();
    data_dim[1] = data_dim[1] + magSf.dimensions()[1];
    T surface_int = Zero;
    forAll(media_patch, face)
    {
        surface_int += data.boundaryField()[media_label][face] *
                       magSf.boundaryField()[media_label][face];
    }
    if (!local)
    {
        surface_int = returnReduce(surface_int, sumOp<T>());
    }
    dimensioned<T> dim_surface_int("surface_int", data_dim, surface_int);
    return dim_surface_int;
}

dimensionedScalar
MacroscaleCompressible::get_surface_area()
{
    const surfaceScalarField& magSf = mesh.magSf();
    scalar surface_int = 0;
    forAll(media_patch, face)
    {
        surface_int += magSf.boundaryField()[media_label][face];
    }
    if (!local)
    {
        surface_int = returnReduce(surface_int, sumOp<double>());
    }
    dimensionedScalar dim_surface_int("surface_int", dimArea, surface_int);
    return dim_surface_int;
}

dimensionedScalar
MacroscaleCompressible::get_volume()
{
    const scalarField& V = mesh.V();
    auto vol_int = sum(V);
    if (!local)
    {
        vol_int = returnReduce(vol_int, sumOp<double>());
    }
    dimensionedScalar dim_volume_int("volume_int", dimVolume, vol_int);
    return dim_volume_int;
}


template <typename T>
dimensioned<T>
MacroscaleCompressible::volume_integrate(
    const GeometricField<T, fvPatchField, volMesh> data)
{
    auto dim_out = data.dimensions() * dimVolume;

    T vol_int = sum(mesh.V() * data.field());

    if (!local)
    {
        vol_int = returnReduce(vol_int, sumOp<T>());
    }
    dimensioned<T> dim_vol_int(dim_out, vol_int);
    return dim_vol_int;
}

void
MacroscaleCompressible::vtk_write_domain()
{
    auto points = averaging_region.points()();

    // Write VTK file header
    vtk_ptr() << "# vtk DataFile Version 3.0\n";
    vtk_ptr() << "Cell Centers with Data (Processor " << Pstream::myProcNo()
              << ")\n";
    vtk_ptr() << "ASCII\n";
    vtk_ptr() << "DATASET UNSTRUCTURED_GRID\n";

    // Write the points to the VTK file
    vtk_ptr() << "POINTS " << points.size() << " float\n";
    forAll(points, pointI)
    {
        const point& p = points[pointI];
        vtk_ptr() << p.x() << " " << p.y() << " " << p.z() << "\n";
    }
    vtk_ptr() << "CELLS " << 1 << " " << 9 << "\n";
    vtk_ptr() << points.size(); // Write number of vertices in the cell
    forAll(points, _p)
    {
        vtk_ptr() << " " << _p; // Write index of each vertex
    }
    vtk_ptr() << "\n";
    vtk_ptr() << "CELL_TYPES " << 1 << "\n";
    vtk_ptr() << "12\n";
    // vtk_ptr() << "CELL_DATA 1\n";
    // vtk_ptr() << "SCALARS temperature float 1\n";
    // vtk_ptr() << "LOOKUP_TABLE default\n";
    // vtk_ptr() << 300.0 << "\n";
    // vtk_ptr() << "VECTORS velocity float\n";
    // vtk_ptr() << Pstream::myProcNo() << " 0.0 0.0\n";
}

void
MacroscaleCompressible::vtk_write_data(word name, dimensionedScalar data)
{
    vtk_ptr() << "CELL_DATA 1\n";
    vtk_ptr() << "SCALARS " << name << " float 1\n";
    vtk_ptr() << "LOOKUP_TABLE default\n";
    vtk_ptr() << data.value() << "\n";
}

void
MacroscaleCompressible::vtk_write_data(word name, dimensionedVector data)
{
    vtk_ptr() << "VECTORS " << name << " float \n";
    vtk_ptr() << data.value().x() << " " << data.value().y() << " "
              << data.value().z() << "\n";
}

} // namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //