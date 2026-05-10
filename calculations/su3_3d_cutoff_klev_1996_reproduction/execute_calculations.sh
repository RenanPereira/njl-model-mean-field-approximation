#!/bin/bash


echo "Script that calculates the integrated cross sections of the SU3 NJL model using the Klevansky \
recipe for the 2-line fermion line integral (B0). The goal of this calculation is to reproduce the results \
presented in the paper:
P. Rehberg, S. P. Klevansky, and J. Hüfner, \
*Elastic scattering and transport coefficients for a quark plasma in SU_f(3) at finite temperatures*, \
Nuclear Physics A 608 (1996) 356-388."
echo ""


# Make neccessary modifications to the codebase to use the Klevansky recipe for the B0 function
./modify_codebase.sh


# Go to project root
cd ../../


# Clean previous build and re-build
make clean
make -j$(nproc)
cp bin/nambuJonaLasinioModel.out calculations/su3_3d_cutoff_klev_1996_reproduction/data
make clean


# Restore to the original
(cd calculations/su3_3d_cutoff_klev_1996_reproduction && ./restore_codebase.sh)


# Use executable to run calculations
cd calculations/su3_3d_cutoff_klev_1996_reproduction/data

filename="integratedCrossSections_setA_TMin0p150_TMax0p250.ini"
./nambuJonaLasinioModel.out use-config-file $filename

rm nambuJonaLasinioModel.out

cd ../../

python3 -m su3_3d_cutoff_klev_1996_reproduction.compute_quark_rel_times
