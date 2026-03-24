from common_utils.shear_viscosity import ShearViscosity


path_input_data_folder = "su3_3d_cutoff_quark_relaxation_times/data/"
path_output_data_folder = "su3_3d_cutoff_transport_coefficients/data/"

# set A
parameter_set = "setA"

method = "COMPLETE_COV"
shear_viscosity = ShearViscosity(
    path_input_data_folder + f'RelaxationTimes_{parameter_set}_{method}_CP0.dat', 
    path_output_data_folder + f'ShearViscosity_{parameter_set}_{method}_CP0.dat'
)

method = "KLEVANSKY"
shear_viscosity = ShearViscosity(
    path_input_data_folder + f'RelaxationTimes_{parameter_set}_{method}_CP0.dat', 
    path_output_data_folder + f'ShearViscosity_{parameter_set}_{method}_CP0.dat'
)

method = "ZHUANG"
shear_viscosity = ShearViscosity(
    path_input_data_folder + f'RelaxationTimes_{parameter_set}_{method}_CP0.dat', 
    path_output_data_folder + f'ShearViscosity_{parameter_set}_{method}_CP0.dat'
)
