import spherepacking

# Delete sphere pack files so runs
spherepacking.clean_directory()

# Determine radii distribution 
radii = spherepacking.SphereRadii(
    n = 100,
    distribution='lognormal',
    mean=0.5,
    stdev=1
)

spheres = spherepacking.Spheres(radii.radii)

domain = spherepacking.Domain(
    spheres=spheres,
    porosity=0.3
)

domain.gen_min_cube()

pack = spherepacking.SpherePack(domain,type = 'Spheres')

pack.gen_input()
pack.run_pack()
pack.read_pack()
pack.save_pack("uniform_spheres")