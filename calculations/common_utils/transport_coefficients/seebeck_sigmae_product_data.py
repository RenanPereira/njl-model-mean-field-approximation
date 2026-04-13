import numpy as np


class SeebeckSigmaeProductData:
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
            self.seebeck_sigmae_product = self.data[:, 1]

        except Exception as e:
            raise ValueError(f"Error loading data from {self.filepath}: {e}")

    def get_temperature(self) -> np.ndarray:
        return self.temperature

    def get_seebeck_sigmae_product(self) -> np.ndarray:
        return self.seebeck_sigmae_product

    def size(self):
        return len(self.temperature)
