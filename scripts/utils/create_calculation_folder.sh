#!/bin/bash

# This script automates the creation of a calculation folder with 
# NO subfolders inside the calculations folder in the root of the project
# It must be executed from the root project using,e.g.:
# (scripts/utils/./create_calculation_folder.sh su3_3d_cutoff_something)

# Validate the input arguments
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <calculation_folder_name>"
    exit 1
fi

calculation_folder_name="$1"

# Create folders
mkdir calculations/${calculation_folder_name}
mkdir calculations/${calculation_folder_name}/data
mkdir calculations/${calculation_folder_name}/plots
mkdir calculations/${calculation_folder_name}/plotting

# Create empty files
touch calculations/${calculation_folder_name}/__init__.py
touch calculations/${calculation_folder_name}/plotting/__init__.py
touch calculations/${calculation_folder_name}/plotting/build_plots.py
touch calculations/${calculation_folder_name}/data/configuration.ini

# Create script that will be used to build plots
cat << EOF > calculations/${calculation_folder_name}/build_plots.sh
#!/bin/bash

echo "Script that generates plots of..."
echo ""

cd ..

python3 -m ${calculation_folder_name}.plotting.build_plots
EOF

# Create script that will be used to execute calculations
cat << EOF > calculations/${calculation_folder_name}/execute_calculations.sh
#!/bin/bash

echo "Script that calculates..."
echo ""

cd .. && cd ..

make -j$(nproc)

cp bin/nambuJonaLasinioModel.out calculations/${calculation_folder_name}/data

cd calculations/${calculation_folder_name}/data

./nambuJonaLasinioModel.out use-config-file CONFIGURATION.ini

rm nambuJonaLasinioModel.out

cd ..
EOF

# Create script that will be used to execute calculations via local main file
cat << 'EOF' > calculations/${calculation_folder_name}/execute_local_main.sh
#!/bin/bash

# This script must be executed from the root folder using: 
# (cd calculations/CALCULATION_FOLDER_NAME && ./execute_local_main.sh )

main_folder_path="calculations/CALCULATION_FOLDER_NAME/data"
main_filename="main.cpp"
number_of_up_folders_to_makefile="2"

../../scripts/utils/switch_src_main_and_run.sh $main_folder_path $main_filename $number_of_up_folders_to_makefile
EOF
sed -i "s/CALCULATION_FOLDER_NAME/${calculation_folder_name}/g" calculations/${calculation_folder_name}/execute_local_main.sh

# Create README.md file
cat << EOF > calculations/${calculation_folder_name}/README.md
# README

## Overview

This folder contains the necessary scripts, data files, and configuration files to...

### Folder Structure
Describe folder structure.

### Prerequisites
- Ensure that...

## How to Use
Explain how the folder should be used.

## Results

Write here some results, append plots, etc.
EOF
