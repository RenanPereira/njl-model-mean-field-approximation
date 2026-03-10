#!/bin/bash

echo "Script that calculates the Cross Sections for the NJL model using different values of temperature and chemical potential"
echo ""

cd .. && cd ..

make -j$(nproc)

cp bin/nambuJonaLasinioModel.out calculations/su3_3d_cutoff_cross_sections_klevansky/data

cd calculations/su3_3d_cutoff_cross_sections_klevansky/data

./nambuJonaLasinioModel.out use-config-file crossSections_T0p215000_CP0p000000.ini
./nambuJonaLasinioModel.out use-config-file crossSections_T0p250000_CP0p000000.ini
./nambuJonaLasinioModel.out use-config-file crossSections_T0p250000_CP0p100000.ini
./nambuJonaLasinioModel.out use-config-file crossSections_T0p250000_CP0p200000.ini
./nambuJonaLasinioModel.out use-config-file crossSections_T0p250000_CP0p300000.ini
./nambuJonaLasinioModel.out use-config-file crossSections_T0p250000_CP0p400000.ini
./nambuJonaLasinioModel.out use-config-file crossSections_T0p250000_CP0p500000.ini

rm nambuJonaLasinioModel.out

cd ..
