import numpy as np

class FixedChemPotTempData:
    def __init__(self, filepath):
        """
        Initialize the class by loading data from the specified file.

        Parameters:
        filepath (str): Path to the data file.
        """
        self.filepath = filepath

        self._load_data()

    def _load_data(self) -> None:
        """
        Private method to load data from the file.
        """
        try:
            self.data = np.loadtxt(self.filepath, skiprows=1)

            self.temperature = self.data[:, 0]
            self.up_quark_effective_mass = self.data[:, 1]
            self.down_quark_effective_mass = self.data[:, 2]
            self.strange_quark_effective_mass = self.data[:, 3]
            self.up_quark_effective_chemical_potential = self.data[:, 4]
            self.down_quark_effective_chemical_potential = self.data[:, 5]
            self.strange_quark_effective_chemical_potential = self.data[:, 6]
            self.pressure = self.data[:, 7]
            self.energy_density = self.data[:, 8]
            self.entropy_density = self.data[:, 9]

        except Exception as e:
            raise ValueError(f"Error loading data from {self.filepath}: {e}")

    def get_temperature(self) -> np.ndarray:
        return self.temperature

    def get_up_quark_effective_mass(self) -> np.ndarray:
        return self.up_quark_effective_mass
    
    def get_down_quark_effective_mass(self) -> np.ndarray:
        return self.down_quark_effective_mass
    
    def get_strange_quark_effective_mass(self) -> np.ndarray:
        return self.strange_quark_effective_mass
    
    def get_up_quark_effective_chemical_potential(self) -> np.ndarray:
        return self.up_quark_effective_chemical_potential
    
    def get_down_quark_effective_chemical_potential(self) -> np.ndarray:
        return self.down_quark_effective_chemical_potential

    def get_strange_quark_effective_chemical_potential(self) -> np.ndarray:
        return self.strange_quark_effective_chemical_potential
    
    def get_pressure(self) -> np.ndarray:
        return self.pressure
    
    def get_energy_density(self) -> np.ndarray:
        return self.energy_density
    
    def get_entropy_density(self) -> np.ndarray:
        return self.entropy_density

    def size(self):
        return len(self.temperature)
