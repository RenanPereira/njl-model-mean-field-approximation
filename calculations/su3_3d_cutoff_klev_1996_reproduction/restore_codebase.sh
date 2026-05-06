#!/bin/bash


echo "Script that restores the codebase after modification made through script modify_codebase.sh"
echo ""


# Go to project root
cd ../../


filepath="src/njl_model/su3_3d_cutoff/SU3NJL3DCutoffMesonPropagators.cpp"

# Backup of the original file
backup_filepath="${filepath}.bak"

# Ensure the backup file exists in the current directory
if [ ! -f "$backup_filepath" ]; then
    echo "Error: Backup file '${backup_file}' not found in the current directory."
    exit 1
fi

# Restore to the original

rm src/njl_model/n_fermion_line_integrals/two_fermion_line_integral_3d_cutoff_klev_recipe.cpp

rm src/njl_model/n_fermion_line_integrals/two_fermion_line_integral_3d_cutoff_klev_recipe.h

cp "$backup_filepath" "$filepath"
rm "$backup_filepath"
