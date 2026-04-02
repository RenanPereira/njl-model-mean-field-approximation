#!/bin/bash

# This script must be executed from the root folder using: 
# (cd calculations/su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp && ./execute_calculations.sh)

echo "Script that calculates the thermodynamics and in medium quark effective masses for the \
SU3 NJL model for different parameter sets as functions of temperature at zero checmical potential."
echo ""

data_folder="calculations/su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp/data"

# Go to project root
cd ../../../

# Clean previous build and re-build
make clean
make -j$(nproc)

cp bin/nambuJonaLasinioModel.out $data_folder

cd $data_folder

filename="SU3NJL3DCutoffFixedChemPotTempThermo_setA_CP0.ini"
./nambuJonaLasinioModel.out use-config-file $filename

filename="SU3NJL3DCutoffFixedChemPotTempThermo_setB_CP0.ini"
./nambuJonaLasinioModel.out use-config-file $filename

filename="SU3NJL3DCutoffFixedChemPotTempThermo_setC_CP0.ini"
./nambuJonaLasinioModel.out use-config-file $filename

filename="SU3NJL3DCutoffFixedChemPotTempThermo_setA_CPCEP.ini"
./nambuJonaLasinioModel.out use-config-file $filename

rm nambuJonaLasinioModel.out

cd ..
