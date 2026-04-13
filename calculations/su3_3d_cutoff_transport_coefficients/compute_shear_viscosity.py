import os
from common_utils.transport_coefficients.shear_viscosity import ShearViscosity


path_input_data_folder = "su3_3d_cutoff_quark_relaxation_times/data/"
path_output_data_folder = "su3_3d_cutoff_transport_coefficients/data/"

parameter_sets = ["setA"]
int_cross_section_methods = ["COMPLETE_COV", "KLEVANSKY", "ZHUANG"]
suffixes = ["CP0", "CPCEP"]

for parameter_set in parameter_sets:
    for method in int_cross_section_methods:
        for suffix in suffixes:
            
            input_filepath = path_input_data_folder + f'RelaxationTimes_{parameter_set}_{method}_{suffix}.dat'
            output_filepath = path_output_data_folder + f'ShearViscosity_{parameter_set}_{method}_{suffix}.dat'
            
            if not os.path.exists(input_filepath):
                print(f"Skipping missing file: {input_filepath}")
                continue
            
            print(f"Calculating the Shear Viscosity based on the file: {input_filepath}")
            print(f"[{parameter_set} | {method} | {suffix}] Processing...")
            
            ShearViscosity(input_filepath, output_filepath)
