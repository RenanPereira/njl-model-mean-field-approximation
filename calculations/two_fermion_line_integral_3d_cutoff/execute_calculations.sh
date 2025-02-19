#!/bin/bash

echo "Script that builds and executes the C++ code that evaluates the two fermion line integral for various scenarios"
echo ""

cd .. && cd ..

make -j$(nproc)

cp bin/nambuJonaLasinioModel.out calculations/two_fermion_line_integral_3d_cutoff/data

cd calculations/two_fermion_line_integral_3d_cutoff/data

./nambuJonaLasinioModel.out use-config-file calculations_B0_vs_k0.ini
./nambuJonaLasinioModel.out use-config-file calculations_B0_vs_k.ini
./nambuJonaLasinioModel.out use-config-file calculations_B0_vs_k0_diff_T.ini
./nambuJonaLasinioModel.out use-config-file calculations_B0_vs_k_diff_T.ini
./nambuJonaLasinioModel.out use-config-file calculations_B0_vs_k0_diff_mu.ini
./nambuJonaLasinioModel.out use-config-file calculations_B0_vs_k_diff_mu.ini
./nambuJonaLasinioModel.out use-config-file calculations_B0_vs_k0_diff_T_mu.ini
./nambuJonaLasinioModel.out use-config-file calculations_B0_vs_k_diff_T_mu.ini

rm nambuJonaLasinioModel.out

cd ..
