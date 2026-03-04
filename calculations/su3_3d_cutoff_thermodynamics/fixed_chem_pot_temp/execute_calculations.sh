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

filename="thermodynamicsSU3NJL3DCutoffFixedChemPotTempSetA.ini"
./nambuJonaLasinioModel.out use-config-file $filename

filename="thermodynamicsSU3NJL3DCutoffFixedChemPotTempSetB.ini"
./nambuJonaLasinioModel.out use-config-file $filename

filename="thermodynamicsSU3NJL3DCutoffFixedChemPotTempSetC.ini"
./nambuJonaLasinioModel.out use-config-file $filename

rm nambuJonaLasinioModel.out

cd ..
