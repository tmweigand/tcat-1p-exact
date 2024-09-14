import numpy as np

class Subdomains:
    """
    Averaging Regions
    """

    def __init__(self, domain, map, overlap = 0.):
        self.domain = domain
        self.map = map
        self.num_subdomains = np.prod(map)
        self.overlap = overlap
        self.length = self.gen_length()

    def gen_boundary_surfaces(self,vol_coords,patch_names):
        """
        Generate the surface patches to get boundary conditions

        Args:
            domain (_type_): _description_
            num_domains (_type_): _description_
        """
        eps = 1.e-3
        coords = []
        names = []
        num_panels = (
            2 * self.map[0] * self.map[1]
            + 2 * self.map[0] * self.map[2]
            + 2 * self.map[1] * self.map[2]
        )

        for n in range(self.num_subdomains):
            index = np.unravel_index(n, self.map)
            for dim in [0,1,2]:
                if index[dim] == 0:
                    _coords = np.copy(vol_coords[n])
                    _coords[dim*2] = self.domain[dim][0] - eps
                    _coords[dim*2+1] = self.domain[dim][0] + eps
                    coords.append(_coords)

                    _name = 'face_'
                    for i in index:
                        _name += str(i)
                    names.append(_name+'_'+str(dim*2))

                if index[dim] == self.map[dim] - 1:
                    _coords = np.copy(vol_coords[n])
                    _coords[dim*2] = self.domain[dim][1] - eps
                    _coords[dim*2+1] = self.domain[dim][1] + eps
                    coords.append(_coords)

                    _name = 'face_'
                    for i in index:
                        _name += str(i)
                    names.append(_name+'_'+str(dim*2+1))

        assert len(coords) == num_panels
        assert len(names) == num_panels
        return coords,names

    def gen_surfaces(self,vol_coords):
        """
        Generate the surface patches

        Args:
            domain (_type_): _description_
            num_domains (_type_): _description_
        """
        eps = 1.e-2
        coords = []
        names = []

        for n in range(self.num_subdomains):
            index = np.unravel_index(n, self.map)
            for dim in [0,1,2]:
                #if index[dim] == 0:
                _coords = np.copy(vol_coords[n])
                _coords[dim*2] -= eps
                _coords[dim*2+1] = _coords[dim*2] + 2*eps
                coords.append(_coords)

                _name = 'face_'
                for i in index:
                    _name += str(i)
                names.append(_name+'_'+str(dim*2))

                #if index[dim] == self.map[dim] - 1:
                _coords = np.copy(vol_coords[n])
                _coords[dim*2] = _coords[dim*2+1] - eps
                _coords[dim*2+1] = _coords[dim*2+1] + eps
                coords.append(_coords)

                _name = 'face_'
                for i in index:
                    _name += str(i)
                names.append(_name+'_'+str(dim*2+1))

        return coords,names


    def gen_length(self):
        """
        Generate length of subdomains
        """
        length = np.zeros(3)
        for dim in [0, 1, 2]:
            length[dim] = (self.domain[dim][1] - self.domain[dim][0]) / self.map[dim]
        return length

    def gen_volumes(self):
        """
        Generate averaging regions boxes and centroids
        """

        coords = []
        names = []
        vertices = []
        bcs = [ [] for _ in range(6)]
        for n in range(self.num_subdomains):
            index = np.unravel_index(n, self.map)
            _coords = np.zeros([6])
            boundary = [0,0,0,0,0,0]
            for dim, ind in enumerate(index):
                pad = self.overlap * self.length[dim]
                if ind == 0 and ind == self.map[dim] - 1:
                    _coords[dim * 2] = self.domain[dim][0]
                    _coords[dim * 2 + 1] = _coords[dim * 2] + self.length[dim] + pad
                    boundary[dim*2] = 1
                    boundary[dim*2+1] = 1
                elif ind == 0:
                    _coords[dim * 2] = self.domain[dim][0]
                    _coords[dim * 2 + 1] = _coords[dim * 2] + self.length[dim] + pad
                    boundary[dim*2] = 1
                elif ind == self.map[dim] - 1:
                    _coords[dim * 2 + 1] = self.domain[dim][1]
                    _coords[dim * 2] = _coords[dim * 2 + 1] - self.length[dim] - pad
                    boundary[dim*2+1] = 1
                else:
                    _coords[dim * 2] = self.domain[dim][0] + ind * self.length[dim]
                    _coords[dim * 2 + 1] = _coords[dim * 2] + self.length[dim] + pad
            coords.append(_coords)

            _name = 'box_'
            for i in index:
                _name += str(i)

            names.append(_name)

            faces = [
                [0,4,7,3],
                [1,2,6,5],
                [4,5,1,0],
                [7,6,2,3],
                [0,1,2,3],
                [4,5,6,7]
            ]

            for n_b,b in enumerate(boundary):
                if b:
                    _face = [8*n+f for f in faces[n_b]]
                    bcs[n_b].append(_face)

            for k in [4,5]:
                ck = _coords[k]
                for j in [2,3]:
                    cj = _coords[j]
                    iL = [0,1]
                    if j == 3:
                        iL = [1,0]
                    for i in iL:
                        ci = _coords[i]
                        _vertex = [ci,cj,ck]
                        vertices.append(_vertex)

        return coords,vertices,names,bcs

    def save_volumes(self, coords, names):
        """
        Print the averaging regions for topoSetDict
        """
        out_file = open("averagingVolumes.in", "w", encoding="utf-8")

        out_file.write("{\n")
        out_file.write("\t name \t interface;\n")
        out_file.write("\t type \t faceSet;\n")
        out_file.write("\t action \t new;\n")
        out_file.write("\t source \t patchToFace;\n")
        out_file.write("\t patch \tmedia;\n")
        out_file.write("}\n\n")

        for n, (c,name) in enumerate(zip(coords,names)):
            index = np.unravel_index(n, self.map)
            out_file.write("{\n")
            out_file.write(f"\t name \t {name};\n")
            out_file.write("\t type \t cellZoneSet;\n")
            out_file.write("\t action \t new;\n")
            out_file.write("\t source \t boxToCell;\n")
            out_file.write("\t boxes\n\t(\n")
            out_file.write(f"\t\t ({c[0]} {c[2]} {c[4]}) ({c[1]} {c[3]} {c[5]})\n")
            out_file.write("\t);\n")
            out_file.write("}\n\n")

            # out_file.write("{\n")
            # out_file.write(f"\t name \t {name};\n")
            # out_file.write("\t type \t cellZoneSet;\n")
            # out_file.write("\t action \t new;\n")
            # out_file.write("\t source \t setToCellZone;\n")
            # out_file.write("\t sourceInfo\n\t{\n")
            # out_file.write(f"\t\t set \t {name+"_set"}; \n")
            # out_file.write("\t}\n")
            # out_file.write("}\n\n")

            out_file.write("{\n")
            out_file.write(f"\t name \t {name+"_ws_set"};\n")
            out_file.write("\t type \t faceSet;\n")
            out_file.write("\t action \t new;\n")
            out_file.write("\t source \t boxToFace;\n")
            out_file.write("\t boxes\n\t(\n")
            out_file.write(f"\t\t ({c[0]} {c[2]} {c[4]}) ({c[1]} {c[3]} {c[5]})\n")
            out_file.write("\t);\n")
            out_file.write("}\n\n")

            out_file.write("{\n")
            out_file.write(f"\t name \t {name+"_ws_set"};\n")
            out_file.write("\t type \t faceSet;\n")
            out_file.write("\t action \t subset;\n")
            out_file.write("\t source \t faceToFace;\n")
            out_file.write("\t set interface;\n")
            out_file.write("}\n\n")

            out_file.write("{\n")
            out_file.write(f"\t name \t {name+"_ws"};\n")
            out_file.write("\t type \t faceZoneSet;\n")
            out_file.write("\t action \t new;\n")
            out_file.write("\t source \t setToFaceZone;\n")
            out_file.write(f"\t faceSet \t {name+"_ws_set"};\n")
            out_file.write("}\n\n")


        out_file.close()

    def save_area(self, coords,names):
        """
        Print the averaging surfaces for topoSetDict
        """
        out_file = open("averagingSurfaces.in", "w", encoding="utf-8")

        for n, (c, name) in enumerate(zip(coords,names)):
            # index = np.unravel_index(n, self.map)
            # name = "face" + str(n)#str(index[0]) + str(index[1]) + str(index[2])
            out_file.write("{\n")
            out_file.write(f"\t name \t {name};\n")
            out_file.write("\t type \t faceZoneSet;\n")
            out_file.write("\t action \t new;\n")
            out_file.write("\t source \t boxToFace;\n")
            out_file.write("\t boxes\n(\n")
            out_file.write(f"\t\t ({c[0]} {c[2]} {c[4]}) ({c[1]} {c[3]} {c[5]})\n")
            out_file.write("\t);\n")
            out_file.write("}\n\n")

        out_file.close()


    def gen_centroids(self, coords):
        """
        Generate the centroids of the averaging regions
        """
        centroid = []
        for c in coords:
            _centroid = (
                c[0] + (c[1] - c[0]) / 2.0,
                c[2] + (c[3] - c[2]) / 2.0,
                c[4] + (c[5] - c[4]) / 2.0,
            )
            centroid.append(_centroid)
        return centroid


    def gen_blockmesh(self,vertices,names,bcs):
        """
        Generate an input script with subdomains for a blockMesh
        """
        out_file = open("blockMeshDict.in", "w", encoding="utf-8")

        out_file.write("vertices\n(\n")
        for v in vertices:
            out_file.write(f"\t({v[0]} {v[1]} {v[2]})\n")
        
        out_file.write(");\n\n")

        out_file.write("blocks\n(\n")
        for n in range(self.num_subdomains):
            out_file.write(f"\t hex (")
            _vertices = np.arange(n*8,(n+1)*8)
            for v in _vertices:
                out_file.write(f"{v} ")
            out_file.write(f") {names[n]} (20 20 20) simpleGrading (1 1 1)\n")

        out_file.write(");\n\n")

        out_file.write("boundary\n(\n")

        bc_names = ['inlet','outlet','top','bottom','front','back']
        bc_types = ['patch','patch','symmetry','symmetry','symmetry','symmetry']


        for n,(bc_name,bc_type) in enumerate(zip(bc_names,bc_types)):
            out_file.write(f"\t {bc_name} \n")
            out_file.write("\t { \n")
            out_file.write(f"\t \t type {bc_type};\n")
            out_file.write("\t \t faces\n \t \t (\n")
            for b in bcs[n]:
                out_file.write(f"\t \t \t ({b[0]} {b[1]} {b[2]} {b[3]})\n")
            out_file.write("\t \t );\n \t }\n \n")

        out_file.write(");\n")

        out_file.close()
