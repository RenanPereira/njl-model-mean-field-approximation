from common_utils.electrical_conductivity import ElectricalConductivity


path_input_data_folder = "su3_3d_cutoff_quark_relaxation_times/data/"
path_output_data_folder = "su3_3d_cutoff_transport_coefficients/data/"

# set A, zero chemical potential
parameter_set = "setA"

method = "COMPLETE_COV"
electrical_conductivity = ElectricalConductivity(
    path_input_data_folder + f'RelaxationTimes_{parameter_set}_{method}_CP0.dat', 
    path_output_data_folder + f'ElectricalConductivity_{parameter_set}_{method}_CP0.dat'
)

method = "KLEVANSKY"
electrical_conductivity = ElectricalConductivity(
    path_input_data_folder + f'RelaxationTimes_{parameter_set}_{method}_CP0.dat', 
    path_output_data_folder + f'ElectricalConductivity_{parameter_set}_{method}_CP0.dat'
)

method = "ZHUANG"
electrical_conductivity = ElectricalConductivity(
    path_input_data_folder + f'RelaxationTimes_{parameter_set}_{method}_CP0.dat', 
    path_output_data_folder + f'ElectricalConductivity_{parameter_set}_{method}_CP0.dat'
)
