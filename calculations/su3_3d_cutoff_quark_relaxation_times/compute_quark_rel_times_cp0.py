from common_utils.quark_relaxation_times import QuarkRelaxationTimes


path_input_data_folder = "su3_3d_cutoff_int_cross_sections/zero_chem_pot/data/"
path_output_data_folder = "su3_3d_cutoff_quark_relaxation_times/data/"

physical_scenario = "zero_chemical_potential_isospin_symmetric"

# set A, zero chemical potential
parameter_set = "setA"

method = "COMPLETE_COV"
quark_relaxation_times = QuarkRelaxationTimes(
    path_input_data_folder, 
    parameter_set, 
    method, 
    physical_scenario,
    path_output_data_folder + f'RelaxationTimes_{parameter_set}_{method}_CP0.dat'
)

method = "KLEVANSKY"
quark_relaxation_times = QuarkRelaxationTimes(
    path_input_data_folder, 
    parameter_set, 
    method, 
    physical_scenario,
    path_output_data_folder + f'RelaxationTimes_{parameter_set}_{method}_CP0.dat'
)

method = "ZHUANG"
quark_relaxation_times = QuarkRelaxationTimes(
    path_input_data_folder, 
    parameter_set, 
    method, 
    physical_scenario,
    path_output_data_folder + f'RelaxationTimes_{parameter_set}_{method}_CP0.dat'
)

# set B, zero chemical potential
parameter_set = "setB"

method = "COMPLETE_COV"
quark_relaxation_times = QuarkRelaxationTimes(
    path_input_data_folder, 
    parameter_set, 
    method, 
    physical_scenario,
    path_output_data_folder + f'RelaxationTimes_{parameter_set}_{method}_CP0.dat'
)

# set C, zero chemical potential
parameter_set = "setC"

method = "COMPLETE_COV"
quark_relaxation_times = QuarkRelaxationTimes(
    path_input_data_folder, 
    parameter_set, 
    method, 
    physical_scenario,
    path_output_data_folder + f'RelaxationTimes_{parameter_set}_{method}_CP0.dat'
)
