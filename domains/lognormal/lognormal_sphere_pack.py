import spherepacking

run_folder = "lognormal"

# Delete sphere pack files so runs
spherepacking.remove_dir(run_folder)

# Determine radii distribution 
radii = spherepacking.SphereRadii(
    n = 100,
    distribution='lognormal',
    mean=0.4,
    stdev=0.2,
    run_folder=run_folder
)

spheres = spherepacking.Spheres(radii.radii)

domain = spherepacking.Domain(
    spheres=spheres,
    porosity=0.3
)

domain.gen_min_cube()

packIO = spherepacking.SpherePackIO(domain,media_type = 'Spheres',run_folder = run_folder, out_folder="data_out")

sp = packIO.generate_sphere_pack(seed = 12, periodic=True)
# packIO.save_pack_txt(sp,"lognormal")
# packIO.save_pack_stl(sp,"uniform")
packIO.save_openfoam(sp,'lognormal_openfoam')