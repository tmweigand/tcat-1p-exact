/*---------------------------------------------------------------------------*\
                Class tcat Declaration
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "objectRegistry.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "volFieldsFwd.H"


#include <iostream>

namespace Foam {

class MacroscaleCompressible {

public:
  const volVectorField &U;
  const volScalarField &p;
  const volScalarField &p_rgh;
  const volScalarField &rho;
  const surfaceScalarField &phi;
  const volScalarField &chem_potential;
  const volScalarField &grav_potential;
  const meshObjects::gravity& g;
  const volScalarField &gh;
  const surfaceScalarField &ghf;
  const dimensionedScalar v_in;
  const dimensionedScalar rho0;
  const dimensionedScalar beta;
  const dimensionedScalar domain_volume;
  const dimensionedScalar mu;
  autoPtr<incompressible::turbulenceModel> turbulence;
  dimensionedScalar e_w;
  dimensionedScalar sauder_mean;

  MacroscaleCompressible(word file_out,const volVectorField &_U, const volScalarField &_p,
                         const volScalarField &_p_rgh,
                         const volScalarField &_rho,
                         const surfaceScalarField &_phi,
                         const volScalarField &_chem_potential,
                         const volScalarField &_grav_potential,
                         const meshObjects::gravity& _g,
                         const volScalarField &_gh, 
                         const surfaceScalarField &_ghf,
                         dimensionedScalar v_in,
                         dimensionedScalar rho0,
                         dimensionedScalar beta,
                         dimensionedScalar domain_volume,
                         dimensionedScalar mu,
                         const int time);

  // Method declaration
  void update(const scalar time);
  dimensionedVector momentum_error();
  dimensionedScalar mass_error();
  dimensionedVector time_int();
  dimensionedVector macro_velocity();
  dimensionedScalar macro_density();
  dimensionedScalar macro_pressure();
  dimensionedScalar macro_chem_potential();
  dimensionedScalar macro_grav_potential();
  dimensionedVector macro_gravity();
  dimensionedTensor macro_stress_tensor(dimensionedScalar macro_pressure, dimensionedVector macro_velocity);
  dimensionedTensor macro_strain_tensor(dimensionedTensor macro_grad_u);
  dimensionedTensor grad_u();
  dimensionedVector grad_p();
  dimensionedVector grad_p_new();
  dimensionedVector grad_e_rho();
  dimensionedVector grad_e_rho_chem();
  dimensionedVector grad_chem_potential();
  dimensionedVector grad_psi();
  dimensionedVector grad_psi_2();
  dimensionedVector div_st();
  dimensionedVector div_vel();
  dimensionedVector surface_int();
  dimensionedScalar get_surface_area();
  dimensionedScalar reynolds(dimensionedScalar sauder_mean, dimensionedScalar rho, dimensionedScalar u);
  tmp<volSymmTensorField> get_stress_tensor();
  tmp<volVectorField> get_div_stress_tensor();

private:
  autoPtr<OFstream> outfile_ptr;
  const fvMesh &mesh;
  label media_label;
  surfaceVectorField normals;
  const polyPatch &media_patch;
  dimensionedVector surface_integrate(const surfaceVectorField &data);
  dimensionedTensor surface_integrate(const surfaceTensorField &data);
};

} // namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
