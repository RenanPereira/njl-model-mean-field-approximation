import numpy as np


class IntegratedCrossSectionData:
    def __init__(self, filepath: str) -> None:
        """
        Initialize the class by loading data from the specified file.

        Parameters:
        filepath (str): Path to the data file.
        """
        self.filepath = filepath
        self.data = np.empty((0, 14))
        
        self.temperature = np.array([])
        self.integrated_cross_section = np.array([])

        self._load_data()

    def _load_data(self) -> None:
        """Private method to load data from the file."""
        try:
            self.data = np.loadtxt(self.filepath, skiprows=1)

            if self.data.shape[1] == 14:
                self.temperature = self.data[:, 0]
                self.integrated_cross_section = self.data[:, 7]
            else:
                raise ValueError(f"Error loading data from {self.filepath}: data file must have 14 columns!")

        except Exception as e:
            raise ValueError(f"Error loading data from {self.filepath}: {e}")

    def get_temperature(self) -> np.ndarray:
        return self.temperature
    
    def get_integrated_cross_section(self) -> np.ndarray:
        return self.integrated_cross_section
