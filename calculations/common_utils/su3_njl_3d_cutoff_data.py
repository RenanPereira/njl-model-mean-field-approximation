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
            self.up_quark_chemical_potential = self.data[:, 4]
            self.down_quark_chemical_potential = self.data[:, 5]
            self.strange_quark_chemical_potential = self.data[:, 6]
            self.pressure = self.data[:, 7]
            self.energy_density = self.data[:, 8]
            self.entropy_density = self.data[:, 9]
            self.up_quark_density = self.data[:, 10]
            self.down_quark_density = self.data[:, 11]
            self.strange_quark_density = self.data[:, 12]

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
    
    def get_up_quark_chemical_potential(self) -> np.ndarray:
        return self.up_quark_chemical_potential
    
    def get_down_quark_chemical_potential(self) -> np.ndarray:
        return self.down_quark_chemical_potential

    def get_strange_quark_chemical_potential(self) -> np.ndarray:
        return self.strange_quark_chemical_potential

    def get_quark_chemical_potential(self, quark_species: str) -> np.ndarray:        
        if quark_species == "up_quark":
            return self.get_up_quark_chemical_potential()
        elif quark_species == "down_quark":
            return self.get_down_quark_chemical_potential()
        elif quark_species == "strange_quark":
            return self.get_strange_quark_chemical_potential()
        else:
            raise ValueError(f"The provided quark_species='{quark_species}' is not valid!")

    def get_baryon_chemical_potential(self) -> np.ndarray:
        chem_pot_u = self.get_up_quark_chemical_potential()
        chem_pot_d = self.get_down_quark_chemical_potential()
        baryon_chem_pot = chem_pot_u + 2.0*chem_pot_d
        return baryon_chem_pot

    def get_pressure(self) -> np.ndarray:
        return self.pressure
    
    def get_energy_density(self) -> np.ndarray:
        return self.energy_density
    
    def get_entropy_density(self) -> np.ndarray:
        return self.entropy_density
    
    def get_up_quark_density(self) -> np.ndarray:
        return self.up_quark_density
    
    def get_down_quark_density(self) -> np.ndarray:
        return self.down_quark_density

    def get_strange_quark_density(self) -> np.ndarray:
        return self.strange_quark_density

    def get_quark_density(self, quark_species: str) -> np.ndarray:        
        if quark_species == "up_quark":
            return self.get_up_quark_density()
        elif quark_species == "down_quark":
            return self.get_down_quark_density()
        elif quark_species == "strange_quark":
            return self.get_strange_quark_density()
        else:
            raise ValueError(f"The provided quark_species='{quark_species}' is not valid!")

    def get_baryon_density(self) -> np.ndarray:
        rho_u = self.get_up_quark_density()
        rho_d = self.get_down_quark_density()
        rho_s = self.get_strange_quark_density()
        baryon_dens = (1.0/3.0)*( rho_u + rho_d + rho_s )
        return baryon_dens
    
    def size(self) -> int:
        return len(self.temperature)
