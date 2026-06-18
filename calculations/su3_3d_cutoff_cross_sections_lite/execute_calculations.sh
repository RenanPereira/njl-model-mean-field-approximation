#!/bin/bash

echo "Script that calculates the Cross Sections for the NJL model."
echo "This calculation is for testing purposes only: a small sample of points are calculated."
echo ""

cd ../.. || exit

make -j"$threads"

cp bin/nambuJonaLasinioModel.out calculations/su3_3d_cutoff_cross_sections_lite/data

cd calculations/su3_3d_cutoff_cross_sections_lite/data

./nambuJonaLasinioModel.out use-config-file crossSections_T0p215000_CP0p000000.ini

rm nambuJonaLasinioModel.out

cd ..

