import numpy as np

class SimulationBox:
    """
    Class to instantiate the bounded box within which the simulation occurs.
    We pass in the no. of cells that would build up an overall box in terms of
    L x B x H. Therefore to have 1000 molecules, we need the box to have
    at least 10x10x10 cells (symmetry is not relevant).

    We first need the dimensions of the unit cell that will be replicated so that
    the box's global coordinates can be established. In doing so, alongwith the
    dimensions of the cell, intersection bounds can be calculated.

    Input: takes a cell object and spawns cells in accordance to the LxBxH stipulation.
    """
    pass

class UnitCell:
    """
    Given an input XYZ coordinate file, will generate a unit cell by calculating
    the minimum cube volume required to avoid overlap with other cells.
    
    For water, its rotation about the Z-axis about the oxygen atom will form a cone.
    This cone can additionally create an approx. spherical volume around it by free
    rotation about its tip. The resulting cube should contain this sphere.

    This object has a unit information attribute that does not change, but is only read
    in order to generate new cells that preserve important properties while still
    introducing random spatial orientations without creating molecular overlaps.
    """
    
    pass

class Molecule:
    # Pass the file path of the XYZ coordinates
    def __init__(self, coordinates_file):
        
        # Read the XYZ coordinates file and create an array of strings, where each
        # element corresponds to a new line in the file. 
        with open(coordinates_file) as file:
            data = file.readlines()
        
        # First line of the XYZ file is the number of atoms. This limits our for loop range.
        n_atoms = int(data[0]) # String -> Int

        atom_coordinates = np.zeros(n_atoms)

        for atom in range(2, n_atoms):
            atom_coordinates