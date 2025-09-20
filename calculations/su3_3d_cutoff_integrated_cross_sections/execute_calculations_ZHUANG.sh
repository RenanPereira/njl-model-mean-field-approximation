#!/bin/bash

echo "Script that calculates the Integrated Cross Sections for the NJL model for different parameter sets \
as functions of temperature and chemical potential. Different methods are used to evaluate the integration of the cross sections."
echo ""

cd .. && cd ..

make clean
make -j$(nproc)

cp bin/nambuJonaLasinioModel.out calculations/su3_3d_cutoff_integrated_cross_sections/data

cd calculations/su3_3d_cutoff_integrated_cross_sections/data

# SetA
filename="integratedCrossSections_setA_TMin0p120000_TMax0p300000.ini"
sed -i 's/approximationMethod = COMPLETE_COV/approximationMethod = ZHUANG/g' $filename
./nambuJonaLasinioModel.out use-config-file $filename
sed -i 's/approximationMethod = ZHUANG/approximationMethod = COMPLETE_COV/g' $filename

rm nambuJonaLasinioModel.out

cd ..
