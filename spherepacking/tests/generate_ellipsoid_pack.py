import spherepacking

# Delete sphere pack files so runs
spherepacking.clean_directory()

# Determine radii distribution 
radii = spherepacking.SphereRadii(
    n = 100,
    distribution='lognormal',
    mean=0.4,
    stdev=0
)

spheres = spherepacking.Spheres(radii.radii)

dim = 2 
factor = 1.5

domain = spherepacking.Domain(
    spheres=spheres,
    porosity=0.3
)

domain.gen_non_cube(dim = dim,factor = factor)

pack = spherepacking.SpherePack(domain,type = "Ellipsoids",dim = dim, factor = factor)

pack.gen_input()
pack.run_pack(output = True)
pack.read_pack()
pack.save_pack()