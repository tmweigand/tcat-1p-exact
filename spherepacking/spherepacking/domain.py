import numpy as np

class Domain:
    """
    Specifications for the Domain
    """

    def __init__(self, spheres, porosity) -> None:
        self.spheres = spheres
        self.porosity = porosity
        self.length = np.zeros(3)

    def gen_min_cube(self, eps=0.0):
        """
        Generate the minimum size cube to fit all spheres
        """
        vol_domain = -self.spheres.volume / (self.porosity - 1.0)
        self.length[:] = np.cbrt(vol_domain) + eps

    def gen_non_cube(self,dim,factor):
        """
        Generate a cube that is longer in a single dimensions
        """
        self.gen_min_cube()
        self.length[dim] = self.length[dim]*factor


    def print_stats(self):
        """
        Print the domain statistics
        """
        print(f"Domain Length: {self.length}")