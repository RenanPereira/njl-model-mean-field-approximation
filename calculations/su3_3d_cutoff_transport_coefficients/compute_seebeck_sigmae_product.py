import os
from common_utils.seebeck_sigmae_product import SeebeckSigmaeProduct


path_quark_rel_times_data_folder = "su3_3d_cutoff_quark_relaxation_times/data/"
path_output_data_folder = "su3_3d_cutoff_transport_coefficients/data/"

parameter_sets = ["setA"]
int_cross_section_methods = ["COMPLETE_COV", "KLEVANSKY", "ZHUANG"]
suffixes = ["CPCEP"]

thermodynamics_filepath = "su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp/data/SU3NJL3DCutoffFixedChemPotTemp_SetA_TMin0p0_TMax0p5_CPU0p318434.dat"

for parameter_set in parameter_sets:
    for method in int_cross_section_methods:
        for suffix in suffixes:
            
            quark_rel_times_filepath = path_quark_rel_times_data_folder + f'RelaxationTimes_{parameter_set}_{method}_{suffix}.dat'
            output_filepath = path_output_data_folder + f'SeebeckSigmaeProduct_{parameter_set}_{method}_{suffix}.dat'
            
            if not os.path.exists(quark_rel_times_filepath):
                print(f"Missing file: {quark_rel_times_filepath}")
            
            print(f"Calculating the SeebeckSigmaeProduct using the files:")
            print(quark_rel_times_filepath)
            print(thermodynamics_filepath)
            print(f"[{parameter_set} | {method} | {suffix}] Processing...")
            
            SeebeckSigmaeProduct(quark_rel_times_filepath, thermodynamics_filepath, output_filepath)
