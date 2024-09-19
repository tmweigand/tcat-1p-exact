import spherepacking

run_folder = "ellipsoid"

# Delete sphere pack files so runs
spherepacking.remove_dir(run_folder)

# Determine radii distribution 
radii = spherepacking.SphereRadii(
    n = 50,
    distribution='normal',
    mean=1,
    stdev=0,
    run_folder=run_folder,
    media_type='ellipsoids'
)

spheres = spherepacking.Ellipsoids(radii.radii)

dim = 2 
factor = 2

domain = spherepacking.Domain(
    spheres=spheres,
    porosity=0.3
)

domain.gen_non_cube(dim = dim,factor = factor)

packIO = spherepacking.SpherePackIO(domain,media_type = 'Ellipsoids',run_folder = run_folder, out_folder="data_out",dim = dim,factor = factor)

sp = packIO.generate_sphere_pack(seed = 12, periodic=True)
# # packIO.save_pack_txt(sp,"uniform")
# # packIO.save_pack_stl(sp,"uniform")
packIO.save_openfoam(sp,'ellipsoid')