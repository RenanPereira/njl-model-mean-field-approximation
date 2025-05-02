import numpy as np


class FirstOrderLineData:
    def __init__(self, filepath):
        """
        Initialize the class by loading data from the specified file.

        Parameters:
        filepath (str): Path to the data file.
        """
        self.filepath = filepath
        self.data = None
        
        self.quark_chem_pot = None
        self.temperature = None
        self.baryon_density_broken = None
        self.baryon_density_restored = None

        self._load_data()

    def _load_data(self):
        """Private method to load data from the file."""
        try:
            self.data = np.loadtxt(self.filepath, skiprows=1)

            self.quark_chem_pot = self.data[:, 0]
            self.temperature = self.data[:, 2]
            self.baryon_density_broken = self.data[:, 4]
            self.baryon_density_restored = self.data[:, 16]

        except Exception as e:
            raise ValueError(f"Error loading data from {self.filepath}: {e}")

    def get_quark_chem_pot(self):
        return self.quark_chem_pot
    
    def get_temperature(self):
        return self.temperature
    
    def get_baryon_density_broken(self):
        return self.baryon_density_broken
    
    def get_baryon_density_restored(self):
        return self.baryon_density_restored
