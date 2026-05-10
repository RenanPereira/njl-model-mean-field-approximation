from common_utils.quark_relaxation_times import QuarkRelaxationTimes


path_input_data_folder = "su3_3d_cutoff_klev_1996_reproduction/data/"
path_output_data_folder = "su3_3d_cutoff_klev_1996_reproduction/data/"


quark_relaxation_times = QuarkRelaxationTimes(
    path_input_data_folder, 
    "setA", 
    "KLEVANSKY", 
    "zero_chemical_potential_isospin_symmetric",
    path_output_data_folder + f'RelaxationTimes_setA_KLEVANSKY_CP0_B0KlevanskyRecipe.dat'
)
