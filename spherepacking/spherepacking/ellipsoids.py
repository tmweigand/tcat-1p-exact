import numpy as np

class Ellipsoids:
    """
    Spheres
    """

    def __init__(self, radii, x=None) -> None:
        self.x = x
        self.radii = radii
        self.volume = 0
        self.gen_volume()

    def gen_volume(self):
        """
        Generate the volume of the spheres
        """
        for r in self.radii:
            self.volume += 4.0 / 3.0 * np.pi * r[0] * r[1] * r[2]

    def print_stats(self):
        """
        Print stats of ellipsoids
        """
        print(f"Volume: {self.volume}")