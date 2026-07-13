from common_utils.quark_relaxation_times import QuarkRelaxationTimes


path_input_data_folder = "su3_3d_cutoff_int_cross_sections/cep_chem_pot/data/"
path_output_data_folder = "su3_3d_cutoff_quark_relaxation_times/data/"

physical_scenario = "finite_chemical_potential_isospin_symmetric"

# set A, CEP chemical potential
parameter_set = "setA"

method = "COMPLETE_COV"
quark_relaxation_times = QuarkRelaxationTimes(
    path_input_data_folder, 
    parameter_set, 
    method, 
    physical_scenario,
    path_output_data_folder + f'RelaxationTimes_{parameter_set}_{method}_CPCEP.dat'
)

# set B, CEP chemical potential
parameter_set = "setB"

method = "COMPLETE_COV"
quark_relaxation_times = QuarkRelaxationTimes(
    path_input_data_folder, 
    parameter_set, 
    method, 
    physical_scenario,
    path_output_data_folder + f'RelaxationTimes_{parameter_set}_{method}_CPCEP.dat'
)

# set C, CEP chemical potential
parameter_set = "setC"

method = "COMPLETE_COV"
quark_relaxation_times = QuarkRelaxationTimes(
    path_input_data_folder, 
    parameter_set, 
    method, 
    physical_scenario,
    path_output_data_folder + f'RelaxationTimes_{parameter_set}_{method}_CPCEP.dat'
)
