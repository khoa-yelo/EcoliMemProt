import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px

from scipy.spatial import ConvexHull

from Bio import PDB


class Coordinate:
    def __init__(self, coords: list):
        self.coords = np.array(coords)
        self.xs = np.array([coordinate[0] for coordinate in coords])
        self.ys = np.array([coordinate[1] for coordinate in coords])
        self.zs = np.array([coordinate[2] for coordinate in coords])

    def subset_coord(self, bool_index):
        sub_coords = []
        for i, coord in enumerate(self.coords):
            if bool_index[i]:
                sub_coords.append(coord)
        return Coordinate(sub_coords)

    def __add__(self, other):
        if isinstance(other, Coordinate):
            return Coordinate(np.concatenate((self.coords, other.coords), axis=0))
        else:
            raise "Invalid addition of Coordinate object"


class OrientedEntity:
    MEMBRANE_TAG = "HETATM"
    AMINO_ACID_TAG = "ATOM"

    def __init__(self, name: str, pdb_file: str):
        self.pdb_file = pdb_file
        self.name = name
        self.entity_coordinates = self.calculate_entity_coordinates()
        self.membrane_coordinates = self.calculate_membrane_coordinates()
        self.coordinates = self.entity_coordinates + self.membrane_coordinates
        self.membrane_num = len(set(self.membrane_coordinates.zs)) / 2

    def calculate_entity_coordinates(self):
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure(self.name, self.pdb_file)
        coordinates = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_id()[0] == " ":
                        for atom in residue:
                            x, y, z = atom.get_coord()
                            coordinates.append([x, y, z])
        return Coordinate(coordinates)

    def calculate_membrane_coordinates(self):
        coordinates = []
        with open(self.pdb_file, "r") as f:
            pdb_text = f.readlines()
            for line in pdb_text:
                if self.MEMBRANE_TAG in line and "DUM" in line:
                    line_vals = tuple(line.strip("\n").split())
                    # fix error where coords are like 100.000-111.000 without space between the 2 numbers
                    tentative_coords = list(line_vals[-3:])
                    cleaned_coords = []
                    for point in reversed(tentative_coords):
                        if len(point) < 9:
                            cleaned_coords.insert(0, point)
                            if len(cleaned_coords) == 3:
                                break
                        else:
                            points = point.split("-")
                            points = list(filter(None, points))
                            for p in reversed(points):
                                if f"-{p}" in point:
                                    p = f"-{p}"
                                cleaned_coords.insert(0, p)
                                if len(cleaned_coords) == 3:
                                    break
                        if len(cleaned_coords) == 3:
                            break
                    line_vals_float = [float(val) for val in cleaned_coords]
                    coordinates.append(line_vals_float)
        return Coordinate(coordinates)

    @property
    def is_transmembrane(self):
        protein_z = self.entity_coordinates.zs
        membrane_z = self.membrane_coordinates.zs

        if min(protein_z) <= min(membrane_z) and max(protein_z) >= max(membrane_z):
            return True
        else:
            return False

    @property
    def is_inward(self):
        protein_z = self.entity_coordinates.zs
        membrane_z = self.membrane_coordinates.zs
        if min(protein_z) <= min(membrane_z) and max(protein_z) < max(membrane_z):
            return True
        else:
            return False

    @property
    def is_outward(self):
        protein_z = self.entity_coordinates.zs
        membrane_z = self.membrane_coordinates.zs
        if min(protein_z) > min(membrane_z) and max(protein_z) > max(membrane_z):
            return True
        else:
            return False

    @property
    def sub_sections(self):
        protein_z = self.entity_coordinates.zs
        membrane_z = self.membrane_coordinates.zs
        inward_index = self.entity_coordinates.zs < min(self.membrane_coordinates.zs)
        outward_index = self.entity_coordinates.zs > max(self.membrane_coordinates.zs)
        trans_index = ~np.logical_or(inward_index, outward_index)
        inward_coord = self.entity_coordinates.subset_coord(inward_index)
        outward_coord = self.entity_coordinates.subset_coord(outward_index)
        trans_coord = self.entity_coordinates.subset_coord(trans_index)

        return {"inward": inward_coord, "trans": trans_coord, "outward": outward_coord}

    @property
    def areas(self):
        areas = {}
        names = ["periplasmic side", "transmembrane", "cytoplasmic side", "total"]
        for entity, name in tuple(
            zip(
                (
                    self.sub_sections["outward"],
                    self.sub_sections["trans"],
                    self.sub_sections["inward"],
                    self.entity_coordinates,
                ),
                names,
            )
        ):
            if len(entity.coords) < 3:
                areas[name] = 0.0
            else:
                entity_points = entity.coords[:, :2]
                hull = ConvexHull(entity_points)
                areas[name] = round(hull.volume / 100.0, 2)  # area in nm^2

        return areas

    def show_footprint(self):
        names = ["periplasmic side", "transmembrane", "cytoplasmic side", "total"]
        i = 0
        padding = 10
        plt.figure(figsize=(8, 8), constrained_layout=True)
        sections = (
            self.sub_sections["outward"],
            self.sub_sections["trans"],
            self.sub_sections["inward"],
            self.entity_coordinates,
        )
        max_coord = max(abs(entity.coords[:, :2]).max() for entity in sections)
        for entity, name in tuple(zip(sections, names)):
            i += 1
            if len(entity.coords) == 0:
                continue
            entity_points = entity.coords[:, :2]
            hull = ConvexHull(entity_points)
            plot_map = {1: 221, 2: 222, 3: 223, 4: 224}
            plt.subplot(plot_map[i])
            plt.title(
                "{0} - area = {1} nm^2".format(name, self.areas[name], fontsize=16)
            )
            plt.plot(entity_points[:, 0], entity_points[:, 1], ".", markersize=1)
            for simplex in hull.simplices:
                plt.plot(entity_points[simplex, 0], entity_points[simplex, 1], "k-")
            n = int(max_coord) + padding
            plt.xlabel("X-axis (Angstrong)")
            plt.ylabel("Y-axis (Angstrong)")
            plt.xlim(-n, n)
            plt.ylim(-n, n)
        plt.show()

    def plot_3d(self, coords: Coordinate = None, color: str = "blue"):
        if coords == None:
            x, y, z = self.coordinates.xs, self.coordinates.ys, self.coordinates.zs
        else:
            x, y, z = coords.xs, coords.ys, coords.zs
        df = pd.DataFrame({"X": x, "Y": y, "Z": z})
        fig = px.scatter_3d(
            df, x="X", y="Y", z="Z", size_max=0.1, title="Interactive 3D Scatter Plot"
        )
        fig.update_traces(
            marker=dict(size=5, color=color)
        )  # Adjust marker size as needed
        fig.show()

    def __repr__(self):
        return self.name


if __name__ == "__main__":
    # Example of usage
    pdb_file = "../data/opm_selected/7nyv.pdb"
    name = "7nyv"
    entity = OrientedEntity(name, pdb_file)
    print(entity.areas)
    entity.show_footprint()
    entity.plot_3d()
