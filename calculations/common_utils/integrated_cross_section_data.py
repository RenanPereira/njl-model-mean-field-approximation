import numpy as np
import os


class IntegratedCrossSectionData:
    def __init__(self, list_filepath: list) -> None:
        self.list_filepath = list_filepath
        
        self.data = np.empty((0, 14))

        self.float_tolerance = 1e-8
        
        self.temperature = np.array([])
        self.up_quark_number = np.array([])
        self.down_quark_number = np.array([])
        self.strange_quark_number = np.array([])
        self.up_antiquark_number = np.array([])
        self.down_antiquark_number = np.array([])
        self.strange_antiquark_number = np.array([])
        self.integrated_cross_section = np.array([])
        self.up_quark_effective_mass = np.array([])
        self.down_quark_effective_mass = np.array([])
        self.strange_quark_effective_mass = np.array([])
        self.up_quark_effective_chemical_potential = np.array([])
        self.down_quark_effective_chemical_potential = np.array([])
        self.strange_quark_effective_chemical_potential = np.array([])

        self._load_data()
        self._sort_data_based_on_temperature()
        self._set_observables()

    def _load_data(self) -> None:
        data_list = []
        for file in self.list_filepath:
            try:
                data  = np.loadtxt(file, skiprows=1)

                if data.shape[1] != 14:
                    raise ValueError(f"Error loading data from {file}: data file must have 14 columns!")
                
                data_list.append(data)    

            except Exception as e:
                raise ValueError(f"Error loading data from {file}: {e}")

        if len(data_list) == 0:
            return

        self.data = np.vstack(data_list) 
    
    def _sort_data_based_on_temperature(self) -> None:
        # Get sorting indices based on the first column and sort based on those
        sort_indices = np.argsort(self.data[:, 0])
        sorted_data = self.data[sort_indices]
        
        if sorted_data.shape[0] == 0:
            return
        
        # Create clean data containing the first row for later comparison
        clean_data = [sorted_data[0]]
        
        # Keep the first occurrence of each value in column 0 (within a small float tolerance)
        for row in sorted_data[1:]:
            if abs(row[0] - clean_data[-1][0]) > self.float_tolerance:
                clean_data.append(row)

        self.data = np.array(clean_data)
    
    def _set_observables(self) -> None:
        self.temperature = self.data[:, 0]
        self.up_quark_number = self.data[:, 1]
        self.down_quark_number = self.data[:, 2]
        self.strange_quark_number = self.data[:, 3]
        self.up_antiquark_number = self.data[:, 4]
        self.down_antiquark_number = self.data[:, 5]
        self.strange_antiquark_number = self.data[:, 6]
        self.integrated_cross_section = self.data[:, 7]
        self.up_quark_effective_mass = self.data[:, 8]
        self.down_quark_effective_mass = self.data[:, 9]
        self.strange_quark_effective_mass = self.data[:, 10]
        self.up_quark_effective_chemical_potential = self.data[:, 11]
        self.down_quark_effective_chemical_potential = self.data[:, 12]
        self.strange_quark_effective_chemical_potential = self.data[:, 13]

    def get_data(self) -> np.ndarray:
        return self.data

    def get_temperature(self) -> np.ndarray:
        return self.temperature

    def get_up_quark_number(self) -> np.ndarray:
        return self.up_quark_number

    def get_down_quark_number(self) -> np.ndarray:
        return self.down_quark_number

    def get_strange_quark_number(self) -> np.ndarray:
        return self.strange_quark_number

    def get_up_antiquark_number(self) -> np.ndarray:
        return self.up_antiquark_number

    def get_down_antiquark_number(self) -> np.ndarray:
        return self.down_antiquark_number

    def get_strange_antiquark_number(self) -> np.ndarray:
        return self.strange_antiquark_number
    
    def get_integrated_cross_section(self) -> np.ndarray:
        return self.integrated_cross_section
    
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

    @classmethod
    def from_matching_files(
        cls, folder: str, parameter_set: str, process: str, method: str
    ) -> "IntegratedCrossSectionData":

        prefix = f'IntegratedCrossSection_{parameter_set}_{process}_{method}_'

        files = [
            (folder + f) for f in os.listdir(folder)
            if f.startswith(prefix)
        ]

        return cls(files)


def read_integrated_cross_section_data(
    folder: str, set: str, process: str, method: str
) -> IntegratedCrossSectionData:

    prefix = f'IntegratedCrossSection_{set}_{process}_{method}_'

    files = [
        (folder + f) for f in os.listdir(folder)
        if f.startswith(prefix)
    ]

    return IntegratedCrossSectionData(files)



class IntegratedCrossSectionDataFromSingleFile:
    def __init__(self, filepath: str) -> None:
        self.filepath = filepath
        self.data = np.empty((0, 14))
        
        self.temperature = np.array([])
        self.integrated_cross_section = np.array([])

        self._load_data()

    def _load_data(self) -> None:
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
