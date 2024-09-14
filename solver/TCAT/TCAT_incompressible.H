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

class MacroscaleIncompressible {

public:
  const volVectorField &U;
  const volScalarField &p;
  const surfaceScalarField &phi;
  autoPtr<incompressible::turbulenceModel> turbulence;

  MacroscaleIncompressible(const volVectorField &_U, const volScalarField &_p,
                      const surfaceScalarField &_phi, const int time);

  // Method declaration
  void update(const scalar time);
  dimensionedVector momentum_error();
  dimensionedScalar mass_error();
  dimensionedVector time_int();
  dimensionedVector grad_p();
  dimensionedVector div_st();
  dimensionedVector div_vel();
  dimensionedVector surface_int();
  dimensionedScalar surface_area_int();
  tmp<volSymmTensorField> get_stress_tensor();
  tmp<volVectorField> get_div_stress_tensor();

private:
  autoPtr<OFstream> outfile_ptr;
  const fvMesh &mesh;
  label media_label;
  surfaceVectorField normals;
  const polyPatch &media_patch;
  dimensionedVector surface_integrate(const surfaceVectorField &data);
};

} // namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
