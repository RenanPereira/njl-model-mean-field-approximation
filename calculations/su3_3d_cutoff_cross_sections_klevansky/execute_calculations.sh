#!/bin/bash

echo "Script that calculates the Cross Sections for the NJL model using different values of temperature and chemical potential"
echo ""

cd .. && cd ..

make -j$(nproc)

cp bin/nambuJonaLasinioModel.out calculations/su3_3d_cutoff_cross_sections_klevansky/data

cd calculations/su3_3d_cutoff_cross_sections_klevansky/data

./nambuJonaLasinioModel.out use-config-file crossSections_T0.215000_CP0.000000.ini
./nambuJonaLasinioModel.out use-config-file crossSections_T0.250000_CP0.000000.ini
./nambuJonaLasinioModel.out use-config-file crossSections_T0.250000_CP0.100000.ini
./nambuJonaLasinioModel.out use-config-file crossSections_T0.250000_CP0.200000.ini
./nambuJonaLasinioModel.out use-config-file crossSections_T0.250000_CP0.300000.ini
./nambuJonaLasinioModel.out use-config-file crossSections_T0.250000_CP0.400000.ini
./nambuJonaLasinioModel.out use-config-file crossSections_T0.250000_CP0.500000.ini

rm nambuJonaLasinioModel.out

cd ..
