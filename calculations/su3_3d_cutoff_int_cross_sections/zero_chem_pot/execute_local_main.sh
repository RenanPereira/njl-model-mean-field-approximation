# This script must be executed from the root folder using: 
# (cd calculations/su3_3d_cutoff_int_cross_sections && ./execute_local_main.sh) 

#!/bin/bash

main_folder_path="calculations/su3_3d_cutoff_int_cross_sections/zero_chem_pot/data"
number_of_up_folders_to_makefile="3"
../../../scripts/utils/switch_src_main_and_run.sh $main_folder_path $number_of_up_folders_to_makefile
