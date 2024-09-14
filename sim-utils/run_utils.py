import generate_subdomains
import generate_functions

domain = [[0, 19.139486], [0, 19.139486], [0, 19.139486]]
overlap = 0.0  # percentage
map = [1, 1, 1]

patch_names = [
    "inlet",
    "outlet",
    "top",
    "bottom",
    "front",
    "back",
]

av = generate_subdomains.Subdomains(domain, map, overlap)

coords, vertices, vol_names, bcs = av.gen_volumes()
av.save_volumes(coords, vol_names)

# av.gen_blockmesh(vertices,vol_names, bcs)

# for s,sn in zip(coords, vol_names):
#     print(s,sn)

# for v in vertices:
#     print(v)

sa, sa_names = av.gen_surfaces(coords)
# av.save_area(sa, sa_names)

# for s,sn in zip(sa,sa_names):
#     print(s,sn)

### Generate Functions for Volume Averages
vol_file_in = "volumeAverageFunctions.in"
vol_variables = [
    "rho",
    "mass_flux",
    "div_u",
    "div_u_rho",
    "div_stress_term_dev",
    "div_stress_term_dev_macro",
    "ddt_u",
    "ddt_rho",
    "ddt_u_rho",
]

generate_functions.volumeIntegrate(vol_file_in, vol_names, vol_variables)

### Generate Functions for Surface Averages
surf_file_in = "surfaceAverageFunctions.in"
surf_variables = ["ten_dot_norm"]

generate_functions.surfaceIntegrate(
    surf_file_in, vol_names, surf_variables, interface=True
)

### Generate Functions for Boundary conditions
surf_file_in = "boundarySurfaceAverageFunctions.in"
surf_variables = [
    "phi",
    "mass_flux",
    "div_u",
    "div_u_rho",
    "stress_term_dev",
    "stress_term_dev_macro",
]

generate_functions.boundarySurfaceIntegrate(surf_file_in, sa_names, patch_names, surf_variables)

# for s,sn in zip(sa,sa_names):
#     print(s,sn)

# centorids = gen_centroids(coords)
# for c, ce in zip(coords, centorids):
#     print(c, ce)
