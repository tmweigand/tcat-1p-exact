/*---------------------------------------------------------------------------*\
                Class tcat Declaration
\*---------------------------------------------------------------------------*/

#include "Pstream.H"
#include "fvCFD.H"
#include "objectRegistry.H"
#include "volFieldsFwd.H"

#include <iostream>

namespace Foam
{

class MacroscaleCompressible
{
public:
    MacroscaleCompressible(bool local,
                           word file_out,
                           const fvMesh& mesh,
                           const int time,
                           dimensionedScalar _mu);

    // Method declaration
    double residual(const volVectorField& U,
                  const volScalarField& p,
                  const volScalarField& p_rgh,
                  const volScalarField& rho,
                  const surfaceScalarField& phi,
                  const meshObjects::gravity& g,
                  const volScalarField& gh,
                  const surfaceScalarField& ghf);

    void update(const int iter,
                const scalar time,
                const volVectorField& U,
                const volScalarField& p,
                const volScalarField& p_rgh,
                const volScalarField& rho,
                const surfaceScalarField& phi,
                const volScalarField& chem_potential,
                const volScalarField& grav_potential,
                const meshObjects::gravity& g,
                const volScalarField& gh,
                const surfaceScalarField& ghf);
    dimensionedScalar time_integral(const volScalarField& arg1);
    dimensionedVector time_integral(const volScalarField& arg1,
                                    const volVectorField& arg2);
    dimensionedScalar get_volume();
    dimensionedScalar get_surface_area();
    dimensionedScalar reynolds(dimensionedScalar sauder_mean,
                               dimensionedScalar rho,
                               dimensionedScalar u);
    tmp<volSymmTensorField> get_stress_tensor(dimensionedScalar mu,
                                              const volVectorField& U);
    tmp<volVectorField> get_div_stress_tensor(dimensionedScalar mu,
                                              const volVectorField& U);

private:
    // Mesh Parameters
    bool local;
    word file_out;

    const fvMesh& mesh;
    label media_label;
    const polyPatch& media_patch;

    dimensionedScalar mu;

    // Outward normal from the w_phase
    surfaceVectorField n_w;
    boundBox averaging_region;
    dimensionedScalar domain_volume;
    point domain_center;
    dimensionedScalar w_volume;
    dimensionedScalar ws_area;
    dimensionedScalar e_w;
    dimensionedScalar e_ws;
    dimensionedScalar sauder_mean;

    // Output
    autoPtr<OFstream> outfile_ptr;
    autoPtr<OFstream> vtk_ptr;

    template <typename T>
    dimensioned<T> average(GeometricField<T, fvPatchField, volMesh> data,
                           dimensionedScalar volume);

    template <typename T>
    dimensioned<T> average(GeometricField<T, fvPatchField, volMesh> data,
                           volScalarField weight);

    dimensionedVector average(const meshObjects::gravity data,
                              volScalarField weight);

    template <typename T1, typename T2>
    dimensioned<T1> divergence(
        tmp<GeometricField<T1, fvPatchField, volMesh> > div_term,
        tmp<GeometricField<T2, fvsPatchField, surfaceMesh> > surface_term);

    template <typename T>
    dimensioned<T> gradient(GeometricField<T, fvPatchField, volMesh> data);

    template <typename T1, typename T2>
    dimensioned<T1>
    gradient(tmp<GeometricField<T1, fvPatchField, volMesh> > grad,
             tmp<GeometricField<T2, fvsPatchField, surfaceMesh> > surface_term);

    template <typename T>
    dimensioned<T>
    surface_integrate(const GeometricField<T, fvsPatchField, surfaceMesh> data);

    template <typename T>
    dimensioned<T>
    volume_integrate(const GeometricField<T, fvPatchField, volMesh> data);

    void vtk_write_domain();
    void vtk_write_data(word name, dimensionedScalar data);
    void vtk_write_data(word name, dimensionedVector data);
};

} // namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
