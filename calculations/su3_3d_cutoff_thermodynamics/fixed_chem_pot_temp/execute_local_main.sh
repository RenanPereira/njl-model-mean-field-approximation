# This script must be executed from the root folder using: 
# (cd calculations/su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp && ./execute_local_main.sh )

#!/bin/bash

main_folder_path="calculations/su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp/data"
number_of_up_folders_to_makefile="3"
../../../tools/utils/switch_src_main_and_run.sh $main_folder_path $number_of_up_folders_to_makefile
