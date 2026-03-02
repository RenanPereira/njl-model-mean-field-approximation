#!/bin/bash

# This script must be executed from the root folder using: 
# (cd calculations/su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp && ./execute_local_main.sh )

main_folder_path="calculations/su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp/data"
main_filename="main.cpp"
number_of_up_folders_to_makefile="3"

../../../scripts/utils/switch_src_main_and_run.sh $main_folder_path $main_filename $number_of_up_folders_to_makefile
