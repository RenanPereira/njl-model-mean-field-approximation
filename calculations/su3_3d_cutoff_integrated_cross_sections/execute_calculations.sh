#!/bin/bash

echo "Script that calculates the Integrated Cross Sections for the NJL model for different parameter sets \
as functions of temperature and chemical potential. Different methods are used to evaluate the integration of the cross sections."
echo ""

cd .. && cd ..

make clean
make -j$(nproc)

cp bin/nambuJonaLasinioModel.out calculations/su3_3d_cutoff_integrated_cross_sections/data

cd calculations/su3_3d_cutoff_integrated_cross_sections/data

./nambuJonaLasinioModel.out use-config-file integratedCrossSections_setA_TMin0.120000_TMax0.300000.ini

rm nambuJonaLasinioModel.out

cd ..
