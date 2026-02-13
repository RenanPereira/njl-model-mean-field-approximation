from common_utils.quark_relaxation_times import QuarkRelaxationTimes


path_input_data_folder = "su3_3d_cutoff_int_cross_sections/zero_chem_pot/data/"
path_output_data_folder = "su3_3d_cutoff_quark_relaxation_times/data/"


# set A, zero chemical potential

quark_relaxation_times = QuarkRelaxationTimes(
    path_input_data_folder, 
    "setA", 
    "COMPLETE_COV", 
    "zero_chemical_potential_isospin_symmetric",
    path_output_data_folder
)

quark_relaxation_times = QuarkRelaxationTimes(
    path_input_data_folder, 
    "setA", 
    "KLEVANSKY", 
    "zero_chemical_potential_isospin_symmetric",
    path_output_data_folder
)

quark_relaxation_times = QuarkRelaxationTimes(
    path_input_data_folder, 
    "setA", 
    "ZHUANG", 
    "zero_chemical_potential_isospin_symmetric",
    path_output_data_folder
)
