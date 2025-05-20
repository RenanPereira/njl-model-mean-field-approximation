import numpy as np


class CrossSectionData:
    def __init__(self, filepath):
        """
        Initialize the class by loading data from the specified file.

        Parameters:
        filepath (str): Path to the data file.
        """
        self.filepath = filepath
        self.data = None
        
        self.sqrt_center_of_mass_energy = None
        self.cross_section = None

        self._load_data()

    def _load_data(self):
        """Private method to load data from the file."""
        try:
            #self.data = np.loadtxt(self.filepath, skiprows=1)
            self.data = np.loadtxt(self.filepath)

            self.sqrt_center_of_mass_energy = self.data[:, 0]
            self.cross_section = self.data[:, 1]

        except Exception as e:
            raise ValueError(f"Error loading data from {self.filepath}: {e}")

    def get_sqrt_center_of_mass_energy(self):
        return self.sqrt_center_of_mass_energy
    
    def get_cross_section(self):
        return self.cross_section
