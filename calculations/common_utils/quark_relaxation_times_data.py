import numpy as np


class QuarkRelaxationTimesData:
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
            self.rel_time_up_quark = self.data[:, 1]
            self.rel_time_down_quark = self.data[:, 2]
            self.rel_time_strange_quark = self.data[:, 3]
            self.rel_time_up_antiquark = self.data[:, 4]
            self.rel_time_down_antiquark = self.data[:, 5]
            self.rel_time_strange_antiquark = self.data[:, 6]
            self.up_quark_effective_mass = self.data[:, 7]
            self.down_quark_effective_mass = self.data[:, 8]
            self.strange_quark_effective_mass = self.data[:, 9]
            self.up_quark_effective_chemical_potential = self.data[:, 10]
            self.down_quark_effective_chemical_potential = self.data[:, 11]
            self.strange_quark_effective_chemical_potential = self.data[:, 12]

        except Exception as e:
            raise ValueError(f"Error loading data from {self.filepath}: {e}")

    def get_temperature(self) -> np.ndarray:
        return self.temperature

    def get_rel_time_up_quark(self) -> np.ndarray:
        return self.rel_time_up_quark

    def get_rel_time_down_quark(self) -> np.ndarray:
        return self.rel_time_down_quark

    def get_rel_time_strange_quark(self) -> np.ndarray:
        return self.rel_time_strange_quark

    def get_rel_time_up_antiquark(self) -> np.ndarray:
        return self.rel_time_up_quark

    def get_rel_time_down_antiquark(self) -> np.ndarray:
        return self.rel_time_down_quark

    def get_rel_time_strange_antiquark(self) -> np.ndarray:
        return self.rel_time_strange_quark

    def get_rel_time(self, quark_species: str) -> np.ndarray:        
        if quark_species == "up_quark":
            return self.get_rel_time_up_quark()
        elif quark_species == "down_quark":
            return self.get_rel_time_down_quark()
        elif quark_species == "strange_quark":
            return self.get_rel_time_strange_quark()
        elif quark_species == "up_antiquark":
            return self.get_rel_time_up_antiquark()
        elif quark_species == "down_antiquark":
            return self.get_rel_time_down_antiquark()
        elif quark_species == "strange_antiquark":
            return self.get_rel_time_strange_antiquark()
        else:
            raise ValueError(f"The provided quark_species='{quark_species}' is not valid!")

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

    def size(self):
        return len(self.temperature)
