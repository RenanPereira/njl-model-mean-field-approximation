from common_utils.quark_relaxation_times import QuarkRelaxationTimes


physical_scenario = "finite_chemical_potential_isospin_symmetric"
path_output_data_folder = "su3_3d_cutoff_quark_relaxation_times/data/"

# set A, CEP chemical potential
path_input_data_folder = "su3_3d_cutoff_int_cross_sections/data/setA_CP0p318436/"
parameter_set = "setA"
method = "COMPLETE_COV"
quark_relaxation_times = QuarkRelaxationTimes(
    path_input_data_folder, 
    parameter_set, 
    method, 
    physical_scenario,
    path_output_data_folder + f'RelaxationTimes_{parameter_set}_{method}_CP0p318436.dat'
)

# set B, CEP chemical potential
path_input_data_folder = "su3_3d_cutoff_int_cross_sections/data/setB_CP0p231030/"
parameter_set = "setB"
method = "COMPLETE_COV"
quark_relaxation_times = QuarkRelaxationTimes(
    path_input_data_folder, 
    parameter_set, 
    method, 
    physical_scenario,
    path_output_data_folder + f'RelaxationTimes_{parameter_set}_{method}_CP0p231030.dat'
)

# set C, CEP chemical potential
path_input_data_folder = "su3_3d_cutoff_int_cross_sections/data/setC_CP0p164012/"
parameter_set = "setC"
method = "COMPLETE_COV"
quark_relaxation_times = QuarkRelaxationTimes(
    path_input_data_folder, 
    parameter_set, 
    method, 
    physical_scenario,
    path_output_data_folder + f'RelaxationTimes_{parameter_set}_{method}_CP0p164012.dat'
)
