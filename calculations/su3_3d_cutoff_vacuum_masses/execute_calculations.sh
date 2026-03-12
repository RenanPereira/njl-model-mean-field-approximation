#!/bin/bash

echo "Script that builds and executes the NJL model to evaluate the effective quark masses for different parameter sets"
echo ""

cd .. && cd ..

make -j$(nproc)

cp bin/nambuJonaLasinioModel.out calculations/su3_3d_cutoff_vacuum_masses/data

cd calculations/su3_3d_cutoff_vacuum_masses/data

./nambuJonaLasinioModel.out use-config-file vacuumSetA.ini
./nambuJonaLasinioModel.out use-config-file vacuumSetB.ini
./nambuJonaLasinioModel.out use-config-file vacuumSetC.ini

rm nambuJonaLasinioModel.out

cd ..
