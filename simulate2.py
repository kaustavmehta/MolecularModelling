import numpy as np

class Molecule:

    def __init__(self, atom_H1, atom_H2, atom_O, gridIndex):
        
        # Coordinates in the box as to where to spawn the molecule
        self.gridIndex = gridIndex
        gridIndex_x = self.gridIndex[0]
        gridIndex_y = self.gridIndex[1]
        gridIndex_z = self.gridIndex[2]

        self.atom_H1 = atom_H1
        atom_H1_x = self.atom_H1[0] * gridIndex_x
        atom_H1_y = self.atom_H1[1] * gridIndex_y
        atom_H1_z = self.atom_H1[2] * gridIndex_z

        self.atom_H2 = atom_H2
        atom_H2_x = self.atom_H2[0] * gridIndex_x
        atom_H2_y = self.atom_H2[1] * gridIndex_y
        atom_H2_z = self.atom_H2[2] * gridIndex_z

        self.atom_O = atom_O
        atom_O_x = self.atom_O[0] * gridIndex_x
        atom_O_y = self.atom_O[1] * gridIndex_y
        atom_O_z = self.atom_O[2] * gridIndex_z


# Read the XYZ coordinates file and create an array of strings, where each
# element corresponds to a new line in the file. 

coordinates_file = "water.xyz"

with open(coordinates_file) as file:
    data = file.readlines()

coordinates = [i.split() for i in data][2:]

atom_H1 = coordinates[0][1:]
atom_O = coordinates[1][1:]
atom_H2 = coordinates[2][1:]

def Spawn(H1, H2, O, gridLocation):

    H1_x, H1_y, H1_z = float(H1[0]) * gridLocation[0], float(H1[1]) * gridLocation[1], float(H1[2]) * gridLocation[2]
    H2_x, H2_y, H2_z = float(H2[0]) * gridLocation[0], float(H2[1]) * gridLocation[1], float(H2[2]) * gridLocation[2]
    O_x, O_y, O_z = float(O[0]) * gridLocation[0], float(O[1]) * gridLocation[1], float(O[2]) * gridLocation[2]

    return [H1_x, H1_y, H1_z], [H2_x, H2_y, H2_z], [O_x, O_y, O_z]

Spawn(atom_H1, atom_H2, atom_O, [1.2, 1.2, 1.2])
Spawn(atom_H1, atom_H2, atom_O, [1.0, 1.0, 1.0])

