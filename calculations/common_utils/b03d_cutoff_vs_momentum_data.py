import numpy as np


class B03DCutoffVsMomentumData:
    def __init__(self, filepath):
        """
        Initialize the class by loading data from the specified file.

        Parameters:
        filepath (str): Path to the data file.
        """
        self.filepath = filepath
        self.data = None
        self.momentum_to_Lambda_ratio = None
        self.Re_B0 = None
        self.Im_B0 = None
        self._load_data()

    def _load_data(self):
        """Private method to load data from the file."""
        try:
            self.data = np.loadtxt(self.filepath, skiprows=1)
            self.momentum_to_Lambda_ratio = self.data[:, 0]
            self.Re_B0 = self.data[:, 1]
            self.Im_B0 = self.data[:, 2]
        except Exception as e:
            raise ValueError(f"Error loading data from {self.filepath}: {e}")

    def get_momentum_to_Lambda_ratio(self):
        """Return the momentum/Lambda values."""
        return self.momentum_to_Lambda_ratio

    def get_Re_B0(self):
        """Return the Re[B0] values."""
        return self.Re_B0

    def get_Im_B0(self):
        """Return the Im[B0] values."""
        return self.Im_B0

    def summary(self):
        """Print a summary of the loaded data."""
        print(f"Data loaded from: {self.filepath}")
        print(f"Number of data points: {len(self.momentum_to_Lambda_ratio)}")
        print(f"momentum/Lambda range: {self.momentum_to_Lambda_ratio.min()} to {self.momentum_to_Lambda_ratio.max()}")
        print(f"Re[B0] range: {self.Re_B0.min()} to {self.Re_B0.max()}")
        print(f"Im[B0] range: {self.Im_B0.min()} to {self.Im_B0.max()}")
