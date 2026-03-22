from common_utils.quark_relaxation_times import QuarkRelaxationTimes


path_input_data_folder = "su3_3d_cutoff_int_cross_sections/cep_chem_pot/data/"
path_output_data_folder = "su3_3d_cutoff_quark_relaxation_times/data/"


# set A, CEP chemical potential
parameter_set = "setA"
physical_scenario = "finite_chemical_potential_isospin_symmetric"

method = "COMPLETE_COV"
quark_relaxation_times = QuarkRelaxationTimes(
    path_input_data_folder, 
    parameter_set, 
    method, 
    physical_scenario,
    path_output_data_folder + f'RelaxationTimes_{parameter_set}_{method}_CPCEP.dat'
)

method = "KLEVANSKY"
quark_relaxation_times = QuarkRelaxationTimes(
    path_input_data_folder, 
    parameter_set, 
    method, 
    physical_scenario,
    path_output_data_folder + f'RelaxationTimes_{parameter_set}_{method}_CPCEP.dat'
)

method = "ZHUANG"
quark_relaxation_times = QuarkRelaxationTimes(
    path_input_data_folder, 
    parameter_set, 
    method, 
    physical_scenario,
    path_output_data_folder + f'RelaxationTimes_{parameter_set}_{method}_CPCEP.dat'
)
