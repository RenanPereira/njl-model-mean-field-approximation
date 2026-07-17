#!/bin/bash

# Parameter passed to the script
calculation="$1"

# This script must be executed from the root folder using: 
# (cd calculations/su3_3d_cutoff_int_cross_sections && ./execute_local_main.sh calculation) 

if [[ "$calculation" == "setA_CP0p000000" || \
      "$calculation" == "setA_CP0p318436" ]]; then

    main_folder_path="calculations/su3_3d_cutoff_int_cross_sections/data/$calculation"
    main_filename="main.cpp"
    number_of_up_folders_to_makefile="2"

    ../../scripts/utils/switch_src_main_and_run.sh $main_folder_path $main_filename $number_of_up_folders_to_makefile

else
    echo "Error: unknown calculation '$calculation'"
    exit 1
fi
