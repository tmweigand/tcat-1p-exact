import spherepacking

run_folder = "uniform"

# Delete sphere pack files so runs
spherepacking.remove_dir(run_folder)

# Determine radii distribution 
radii = spherepacking.SphereRadii(
    n = 100,
    distribution='normal',
    mean=1,
    stdev=0,
    run_folder=run_folder,
    media_type="spheres"
)

spheres = spherepacking.Spheres(radii.radii)

domain = spherepacking.Domain(
    spheres=spheres,
    porosity=0.3
)

domain.gen_min_cube()

packIO = spherepacking.SpherePackIO(domain,media_type = 'Spheres',run_folder = run_folder, out_folder="data_out")

sp = packIO.generate_sphere_pack(seed = 12, periodic=True)
# packIO.save_pack_txt(sp,"uniform")
# packIO.save_pack_stl(sp,"uniform")
packIO.save_openfoam(sp,'uniform_openfoam')