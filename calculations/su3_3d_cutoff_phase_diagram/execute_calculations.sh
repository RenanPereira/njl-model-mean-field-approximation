#!/bin/bash

echo "Script that builds and executes the NJL model phase diagram for different parameter sets"
echo ""

cd .. && cd ..

make -j$(nproc)

cp bin/nambuJonaLasinioModel.out calculations/su3_3d_cutoff_phase_diagram/data

cd calculations/su3_3d_cutoff_phase_diagram/data

./nambuJonaLasinioModel.out use-config-file firsOrderLineSetA.ini
./nambuJonaLasinioModel.out use-config-file firsOrderLineSetB.ini
./nambuJonaLasinioModel.out use-config-file firsOrderLineSetC.ini

rm nambuJonaLasinioModel.out

cd ..
