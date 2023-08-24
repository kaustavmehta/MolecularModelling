import numpy as np
from scipy.spatial.transform import Rotation as R

class Molecule:

    def __init__(self, coordinates_file=False, initialPos=False) -> None:
        
        if initialPos != False:
            self.initialPos = initialPos # Private data member that is read-only
            self.molPos = self.initialPos # Molecular coordinates are initialized with the initial position.

        # If coordinates_file is not False, we open the file and read its contents
        if coordinates_file != False:
            # Use context manager to safely read the file
            
            with open(coordinates_file) as file:
                data = file.readlines()
                data = [i.split() for i in data]
                # print(data[2:])

            self.nAtoms = int(data[0][0]) # Read the 1st line for the number of atoms
            
            # Array of tuples that stores coordinates as an nparray (1,3) for every element
            coords = [(str(atom[0][0]), np.array([atom[1:]])) for atom in data[2:]]
            # print(coords)
            
            elements, posMat = np.array([coords[0][0]]), np.array(coords[0][1])

            for atoms in coords[1:]:
                # print(atoms)
                element, pos = atoms[0], atoms[1]
                elements = np.vstack([elements, [element]]) # Create a column vector from the element order
                posMat = np.vstack([posMat, pos]) # Create a 3xN position matrix

            self.elements = elements
            self.initialPos = posMat.astype(np.float32)
            self.molPos = self.initialPos
        
    def setPos(self, pos:np.ndarray) -> None:
        """Sets the molecule position using a supplied position matrix (nparray)."""
        self.molPos = pos

    def rotate(self, angleX, angleY, angleZ, degrees=True):
        """Rotates by the given per-coordinate angle (rad) and returns a transformed position matrix."""
        
        # Creates a rotation matrix using scipy.spatial.transform by stacking each of the coordinate rotations
        rotVec = R.from_rotvec([[angleX, 0, 0],
                                [0, angleY, 0],
                                [0, 0, angleZ]], degrees=degrees)
        
        rotVec.as_matrix()

        result = rotVec.apply(self.molPos)
        self.setPos(result)
    
    def randRotate(self, state=123) -> None:
        """Rotates by an angle randomly sampled from a uniform distribution."""

        # Generate a random rotation matrix in three dimensions
        rotVec = R.random(1, random_state=state)
        rotVec.as_matrix()
        result = rotVec.apply(self.molPos) # Compute rotation transformation on the current position matrix
        
        self.setPos(result) # Set new coordinates after random rotation

water = Molecule(coordinates_file='water.xyz')

water.molPos.shape
water.elements
water.randRotate()
water.molPos
water.initialPos