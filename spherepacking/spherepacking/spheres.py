import numpy as np

class Spheres:
    """
    Spheres
    """

    def __init__(self, radii, x=None) -> None:
        self.x = x
        self.radii = radii
        self.n = len(radii)
        self.volume = 0
        self.gen_volume()

    def gen_volume(self):
        """
        Generate the volume of the spheres
        """
        for r in self.radii:
            self.volume += 4.0 / 3.0 * np.pi * r * r * r

    def print_stats(self):
        """
        Print stats of radii
        """
        print(f"Volume: {self.volume}")
        print(f"Mean: {np.mean(self.radii)}")
        print(f"St Dev: {np.std(self.radii)}")