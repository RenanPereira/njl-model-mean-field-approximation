#!/bin/bash

main_folder_path="calculations/su3_3d_cutoff_int_cross_sections/cep_chem_pot/data"
main_filename="main.cpp"
number_of_up_folders_to_makefile="3"

../../../scripts/utils/switch_src_main_and_run.sh $main_folder_path $main_filename $number_of_up_folders_to_makefile
