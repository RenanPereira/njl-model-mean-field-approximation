import numpy as np


class ShearViscosityData:
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
            self.shear_viscosity = self.data[:, 1]

        except Exception as e:
            raise ValueError(f"Error loading data from {self.filepath}: {e}")

    def get_temperature(self) -> np.ndarray:
        return self.temperature

    def get_shear_viscosity(self) -> np.ndarray:
        return self.shear_viscosity

    def size(self):
        return len(self.temperature)
