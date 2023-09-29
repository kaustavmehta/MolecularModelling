import numpy as np
from scipy.spatial.transform import Rotation as R
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform

class Molecule:

    # coordinates_file: 
    #   Default: False
    #   This argument can take a valid file path to a XYZ coordinate file to use those coordinates to instantiate the molecule 
    # initialPos:
    #   Default: False
    #   This argument can take a numpy array of [X, Y, Z] coordinates to instantiate the molecule with these coordinates.
    def __init__(self, coordinates_file=False, initialPos=False) -> None:
        
        # If we supply an initial position for the molecule, update the class attributes to be this initial position
        if initialPos != False:
            self.initialPos = initialPos # Private data member that is read-only
            self.molPos = self.initialPos # Molecular coordinates are initialized with the initial position.

        # If coordinates_file is not False, we open the file and read its contents
        if coordinates_file != False:
            # Use context manager to safely read the file
            with open(coordinates_file) as file:
                data = file.readlines() # Read all the data line by line
                data = [i.split() for i in data] # List comprehension to split by whitespace

            self.nAtoms = int(data[0][0]) # Read the 1st line for the number of atoms
            
            # Array of tuples that stores coordinates as a numpy array of dimension (1,3) for every element
            coords = [(str(atom[0][0]), np.array([atom[1:]])) for atom in data[2:]] # List comprehension
            
            # Array of element names and an array of corresponding coordinates
            elements, posMat = np.array([coords[0][0]]), np.array(coords[0][1])

            for atoms in coords[1:]:
                # Unpack element name and corresponding coordinate
                element, pos = atoms[0], atoms[1]
                # Vertically stack the array to create column vector of element names
                elements = np.vstack([elements, [element]]) # Create a column vector from the element order
                posMat = np.vstack([posMat, pos]) # Create a 3xN position matrix

            self.elements = elements
            self.initialPos = posMat.astype(np.float32) # Can lower the precision to save on memory consumption
            self.molPos = self.initialPos
    
    # Setter method for the molecule's current position
    def setPos(self, pos:np.ndarray) -> None:
        """Sets the molecule position using a supplied position matrix (nparray)."""
        self.molPos = pos

    # Rotation function to rotate by an arbitrary amount in either degree or radians
    def rotate(self, angleX, angleY, angleZ, degrees=True):
        """Rotates by the given per-coordinate angle (rad) and returns a transformed position matrix."""
        
        # Creates a rotation matrix using scipy.spatial.transform by stacking each of the coordinate rotations
        # Rotation matrix is generated using quaternions implicitly
        rotVec = R.from_rotvec([[angleX, 0, 0],
                                [0, angleY, 0],
                                [0, 0, angleZ]], degrees=degrees)
        
        # Get the generated rotation quaternion as a matrix that can be used as a transformation
        rotVec.as_matrix()

        # Apply the rotation matrix to the psition matrix of the molecule to get the final coordinates of the
        # rotated position
        result = rotVec.apply(self.molPos)
        # Use the setter method to update the molecule's position
        self.setPos(result)
    
    def randRotate(self) -> None:
        """Rotates by an angle randomly sampled from a uniform distribution."""

        # Generate a random rotation matrix in three dimensions
        rotVec = R.random(1) # Generate one randomly generated rotational quaternion
        rotVec.as_matrix() # Represent the rotational transformation as a matrix
        result = rotVec.apply(self.molPos) # Compute rotation transformation on the current position matrix
        
        self.setPos(result) # Set new coordinates after random rotation

    # Translate a supplied Molecule object by a random amount in either X,Y,Z directions.
    # Both minimum and maximum bounds can be set which set the limits of the uniform distribution from which the
    # values are sampled from.
    def randTranslate(self, minBound=-1.0, maxBound=1.0, inplace=True) -> None:
        """
        ## Ensuring non-intersection of the next molecule relative to the current molecule/basis:

        We use the O atom as the center to spawn the next molecule. This
        is because we define a bounding spherical volume from R=1.2084 A that contains
        the current water molecule where there is a guarantee that it won't breach
        the bounds of this sphere. Therefore, given any rotation operation on the next molecule,
        we can be sure about its bounds. 
        
        We can then use this radius to construct 2R, the diameter. Let X be the point
        where the next molecule will spawn. Then, as long as the translation vector applied on
        the original position vector of the O atom, X0, yields a vector 
        whose norm is at least (2R)^2, we can be guaranteed that the next molecule will not
        spawn within the intersection radius of the bounding sphere created by the current molecule.
        ----
        ## Ensuring next molecule does not intersect with any molecule in space:
        
        We construct a valid distance matrix based on N atom coordinates. We then construct a
        hypothetical distance matrix based on the N+1 atom coordinates and run a vectorized boolean
        query to check if the distance is strictly greater than 2R for every point. This makes the
        N+1 atom coordinate a valid entry to the system.
        ----
        """
        # Use a default RNG function that is generated when numpy is imported. We convert this to an uniform distribution
        # from which we can sample using the minimum and maximum bounds as the upper and lower limits of this distribution.
        # We also want a column vector of dimension (1,3) because we want to translate every atom in the molecule by the same amount.
        randTranslationVec = np.random.default_rng().uniform(minBound, maxBound, (1,3))
        # Matrix element-wise addition
        result = self.molPos + randTranslationVec

        # Allows for an inplace update of the supplied molecule's position using the randomly generated translation vector
        if inplace == True:
            self.setPos(result)
        # If not inplace=True, we return the coordinates of the randomly translated molecule that the user would have to
        # manually set the position of the molecule using. 
        else:
            return result


class System:

    def __init__(self, ) -> None:
        pass        



water = Molecule(coordinates_file='water.xyz')

water.molPos.shape
water.elements[0][0]
water.randRotate()
water.molPos.astype(np.str_)[0][0]
water.initialPos
water.nAtoms

def MolToXYZ(Molecule):
    positions = Molecule.molPos
    elements = Molecule.elements

    lines = []

    #lines.append(f"{Molecule.nAtoms} \n")
    #lines.append("XYZ file generated by SimulationScriptv1. \n")

    for element, position in zip(elements, positions.astype(np.str_)):
        lines.append(f"{element[0]}    {position[0]}   {position[1]}   {position[2]}\n")

    return lines

# Pass an array of Molecule objects to write the final XYZ file.
def writeXYZ(Molecules):

    nAtoms = len(Molecules) * 3 # Counter to track the number of total atoms (scales with no. of molecules)
    _lines = []

    for molecule in Molecules:
        _lines = _lines + MolToXYZ(molecule)
    
    lines = [f'{nAtoms}\n', 'Molecule auto-generated by script.\n']
    lines = lines + _lines
    
    return lines

water1 = Molecule(coordinates_file='water.xyz')
water1.randRotate()
water1.randTranslate(minBound=np.square(1.21), maxBound=np.square(1.21*100), inplace=True)
water1.molPos

water2 = Molecule(coordinates_file='water.xyz')
water2.randRotate()
water2.randTranslate(minBound=np.square(1.21), maxBound=np.square(1.21*100), inplace=True)
water2.molPos

writeXYZ([water1, water2])

with open('testWaters1.xyz', 'w') as file:
    file.writelines(writeXYZ([water1, water2]))

## Testing 10 molecules

molecules = [Molecule(coordinates_file='water.xyz') for i in range(100000)]

molecules[0].molPos[1].shape

for molecule in molecules:
    molecule.randRotate()
    molecule.randTranslate(minBound=np.square(2.1), maxBound=np.square(1.21*1000), inplace=True)

system_coordinates = np.vstack([molecule.molPos[1] for molecule in molecules])
dist_vector = pdist(system_coordinates)
dist_matrix = squareform(dist_vector)

indices = np.transpose(np.nonzero(dist_matrix < 2.93))


writeXYZ(molecules)

with open('testWaters100000.xyz', 'w') as file:
    file.writelines(writeXYZ(molecules))