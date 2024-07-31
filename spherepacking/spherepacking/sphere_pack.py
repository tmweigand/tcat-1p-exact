import numpy as np
import subprocess
import pygmsh

from .spheres import Spheres
from .ellipsoids import Ellipsoids

class SpherePack:
    """
    Class for running the sphere pack code
    """

    def __init__(self, domain, type, dim = None, factor = None) -> None:
        self.domain = domain
        self.media_type = type
        self.dim = dim
        self.factor = factor
        self.media = None
        self.length = np.zeros(3)
        self.porosity = 0
        self.n_spheres = 0

    def gen_input(self, seed=333, contraction_rate=1.328910e-2):
        """
        Generate 'generation.conf' for the sphere packing code
        """
        out_file = open("generation.conf", "w", encoding="utf-8")
        out_file.write(f"Particles count: {self.domain.spheres.n}\n")
        out_file.write(
            f"Packing size: {self.domain.length[0]} {self.domain.length[1]} {self.domain.length[2]}\n"
        )
        out_file.write("Generation start: 1\n")
        out_file.write(f"Seed: {seed}\n")
        out_file.write("Steps to write: 0\n")
        out_file.write("Boundaries mode: 1\n")
        out_file.write(f"Contraction rate: {contraction_rate} \n")
        out_file.write("Generation mode: 1\n")
        out_file.close()

    def run_pack(self, output=False):
        """
        Run the sphere pack code
        """
        if not output:
            subprocess.run(
                args=["./PackingGeneration.exe", "-fba"],
                check=False,
            )
        else:
            subprocess.run(
                args=["./PackingGeneration.exe", "-fba"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=False,
            )

    def read_pack(self):
        """
        Read the results of the sphere pack
        """

        with open("packing.nfo") as nfo_file:
            head = [next(nfo_file) for x in range(5)]

        self.n_spheres = int(head[0].split(":")[1])

        split = head[1].split(":")[1].split(" ")
        for n, dim in enumerate([1, 2, 3]):
            self.length[n] = float(split[dim])

        theory_poro = float(head[2].split(":")[1])
        act_poro = float(head[3].split(":")[1].split(" ")[1])
        scale_fac = pow((1.0 - act_poro) / (1.0 - theory_poro), 1.0 / 3.0)

        sphere_data = np.fromfile("packing.xyzd")

        x = np.zeros([self.n_spheres,3])
        r = np.zeros(self.n_spheres)

        for n, i in enumerate(range(0, self.n_spheres * 4, 4)):
            x[n,0] = sphere_data[i]
            x[n,1] = sphere_data[i + 1]
            x[n,2] = sphere_data[i + 2]
            r[n] = sphere_data[i + 3] * scale_fac / 2

        if self.media_type == 'Spheres':
            self.media = Spheres(r, x)
            self.porosity = 1 - self.media.volume / np.prod(self.length)
        
        if self.media_type == "Ellipsoids":
            
            radii = np.zeros([self.n_spheres,3])
            for n in range(self.n_spheres):
                for d in [0,1,2]:
                    if d == self.dim:
                        radii[n,d] = r[n]/self.factor
                        x[n,d] = x[n,d]/self.factor
                    else:
                        radii[n,d] = r[n]

            self.media = Ellipsoids(radii,x)
            self.porosity = 1 - self.media.volume / np.prod(self.length)

    def print_stats(self):
        """
        Print the statistics
        """
        print(f"Length {self.length}")
        print(f"Sphere Volume: {self.media.volume}")
        print(f"Domain Volume: {np.prod(self.length)}")
        print(f"Porosity: {self.porosity}")

    def save_domain(self):
        """
        Save the domain as vtk for viz
        """
        min_length = np.min(self.length)
        with pygmsh.geo.Geometry() as geom:
            geom.add_box(
                0.0,
                self.length[0],
                0.0,
                self.length[1],
                0.0,
                self.length[2],
                mesh_size=0.5 * min_length,
            )
            mesh = geom.generate_mesh()
            mesh.write("domain.stl")

    def save_pack(self):
        """
        Save the pack as a stl
        """
        self.save_domain()

        if self.media_type == 'Spheres':
            with pygmsh.geo.Geometry() as geom:
                for n in range(self.n_spheres):
                    geom.add_ball(
                        [self.media.x[n,0], self.media.x[n,1], self.media.x[n,2]],
                        self.media.radii[n],
                        mesh_size=0.5,
                    )
                mesh = geom.generate_mesh()
                mesh.write("spheres.stl")

        elif self.media_type == 'Ellipsoids':
            with pygmsh.geo.Geometry() as geom:
                for n in range(self.n_spheres):
                    geom.add_ellipsoid(
                        [self.media.x[n,0], self.media.x[n,1], self.media.x[n,2]],
                        [
                            self.media.radii[n,0],
                            self.media.radii[n,1],
                            self.media.radii[n,2],
                        ],
                        mesh_size=0.5,
                    )
                mesh = geom.generate_mesh()
                mesh.write("ellipsoids.stl")
