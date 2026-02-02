#!/bin/bash

echo "Script that calculates the Integrated Cross Sections for the NJL model for different parameter sets \
as functions of temperature and chemical potential. Different methods are used to evaluate the integration of the cross sections."
echo ""

# Go to project root
cd .. && cd .. && cd ..

# Clean previous build and re-build
make clean
make -j$(nproc)

cp bin/nambuJonaLasinioModel.out calculations/su3_3d_cutoff_int_cross_sections/zero_chem_pot/data

cd calculations/su3_3d_cutoff_int_cross_sections/zero_chem_pot/data

# SetA
filename="integratedCrossSections_setA_TMin0p120000_TMax0p300000.ini"
sed -i 's/approximationMethod = COMPLETE_COV/approximationMethod = KLEVANSKY/g' $filename
./nambuJonaLasinioModel.out use-config-file $filename
sed -i 's/approximationMethod = KLEVANSKY/approximationMethod = COMPLETE_COV/g' $filename

rm nambuJonaLasinioModel.out

cd ..
