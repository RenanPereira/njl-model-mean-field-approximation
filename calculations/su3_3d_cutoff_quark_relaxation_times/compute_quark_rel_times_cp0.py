from common_utils.quark_relaxation_times import QuarkRelaxationTimes


path_output_data_folder = "su3_3d_cutoff_quark_relaxation_times/data/"
physical_scenario = "zero_chemical_potential_isospin_symmetric"

# set A, zero chemical potential
path_input_data_folder = "su3_3d_cutoff_int_cross_sections/data/setA_CP0p000000/"
parameter_set = "setA"
method = "COMPLETE_COV"
quark_relaxation_times = QuarkRelaxationTimes(
    path_input_data_folder, 
    parameter_set, 
    method, 
    physical_scenario,
    path_output_data_folder + f'RelaxationTimes_{parameter_set}_{method}_CP0.dat'
)

path_input_data_folder = "su3_3d_cutoff_int_cross_sections/data/setA_CP0p000000/"
parameter_set = "setA"
method = "KLEVANSKY"
quark_relaxation_times = QuarkRelaxationTimes(
    path_input_data_folder, 
    parameter_set, 
    method, 
    physical_scenario,
    path_output_data_folder + f'RelaxationTimes_{parameter_set}_{method}_CP0.dat'
)

path_input_data_folder = "su3_3d_cutoff_int_cross_sections/data/setA_CP0p000000/"
parameter_set = "setA"
method = "ZHUANG"
quark_relaxation_times = QuarkRelaxationTimes(
    path_input_data_folder, 
    parameter_set, 
    method, 
    physical_scenario,
    path_output_data_folder + f'RelaxationTimes_{parameter_set}_{method}_CP0.dat'
)

# set B, zero chemical potential
path_input_data_folder = "su3_3d_cutoff_int_cross_sections/data/setB_CP0p000000/"
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
path_input_data_folder = "su3_3d_cutoff_int_cross_sections/data/setC_CP0p000000/"
parameter_set = "setC"
method = "COMPLETE_COV"
quark_relaxation_times = QuarkRelaxationTimes(
    path_input_data_folder, 
    parameter_set, 
    method, 
    physical_scenario,
    path_output_data_folder + f'RelaxationTimes_{parameter_set}_{method}_CP0.dat'
)
