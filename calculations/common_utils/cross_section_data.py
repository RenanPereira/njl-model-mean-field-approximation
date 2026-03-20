import numpy as np


PROCESS_LATEX = {
    "UUUU": r"uu \rightarrow uu",
    "UDUD": r"ud \rightarrow ud",
    "USUS": r"us \rightarrow us",
    "SSSS": r"ss \rightarrow ss",
    "UUBarUUBar": r"u\bar{u} \rightarrow u\bar{u}",
    "UUBarDDBar": r"u\bar{u} \rightarrow d\bar{d}",
    "UUBarSSBar": r"u\bar{u} \rightarrow s\bar{s}",
    "UDBarUDBar": r"u\bar{d} \rightarrow u\bar{d}",
    "USBarUSBar": r"u\bar{s} \rightarrow u\bar{s}",
    "SUBarSUBar": r"s\bar{u} \rightarrow s\bar{u}",
    "SSBarUUBar": r"s\bar{s} \rightarrow u\bar{u}",
    "SSBarSSBar": r"s\bar{s} \rightarrow s\bar{s}",
    "UBarUBarUBarUBar": r"\bar{u}\bar{u} \rightarrow \bar{u}\bar{u}",
    "UBarDBarUBarDBar": r"\bar{u}\bar{d} \rightarrow \bar{u}\bar{d}",
    "UBarSBarUBarSBar": r"\bar{u}\bar{s} \rightarrow \bar{u}\bar{s}",
    "SBarSBarSBarSBar": r"\bar{s}\bar{s} \rightarrow \bar{s}\bar{s}",
}


def process_to_latex(process: str) -> str:
    try:
        return fr'${PROCESS_LATEX[process]}$'
    except KeyError:
        raise ValueError(f"Unknown process '{process}'")


def process_to_ylabel_latex(process: str) -> str:
    try:
        y_label = r"\overline{\sigma}_{" + PROCESS_LATEX[process] + r"} \, [\mathrm{GeV}^{-2}]"
        return fr'${y_label}$'
    except KeyError:
        raise ValueError(f"Unknown process '{process}'")

class CrossSectionData:
    def __init__(self, filepath):
        """
        Initialize the class by loading data from the specified file.

        Parameters:
        filepath (str): Path to the data file.
        """
        self.filepath = filepath
        self.data = np.empty((0, 2))
        
        self.sqrt_center_of_mass_energy = np.array([])
        self.cross_section = np.array([])

        self._load_data()

    def _load_data(self):
        """Private method to load data from the file."""
        try:
            self.data = np.loadtxt(self.filepath, skiprows=1)

            if self.data.ndim == 2:
                self.sqrt_center_of_mass_energy = self.data[:, 0]
                self.cross_section = self.data[:, 1]
            else:
                raise ValueError(f"Error loading data from {self.filepath}: data file must have two columns!")

        except Exception as e:
            raise ValueError(f"Error loading data from {self.filepath}: {e}")

    def get_sqrt_center_of_mass_energy(self):
        return self.sqrt_center_of_mass_energy
    
    def get_cross_section(self):
        return self.cross_section
