#!/bin/bash

echo "Script that calculates the Integrated Cross Sections for the NJL model for different parameter sets \
as functions of temperature and chemical potential. Different methods are used to evaluate the integration of the cross sections."
echo ""

# Go to project root
cd ../../.. || exit

# Change code before compilation: in this case the value of M2 is arbitrary, 
# the switch suggested below will be made permanent after proper testing
filename_one_fermion_line_implementation="src/njl_model/n_fermion_line_integrals/one_fermion_line_integral_3d_cutoff.cpp"
sed -i 's/double M2 = (1-1E-4)\*cutoff;/double M2 = M1;/g' $filename_one_fermion_line_implementation

# Clean previous build and re-build
make clean
make -j$(nproc)

# Undo the change made above after the build
sed -i 's/double M2 = M1;/double M2 = (1-1E-4)\*cutoff;/g' $filename_one_fermion_line_implementation

cp bin/nambuJonaLasinioModel.out calculations/su3_3d_cutoff_int_cross_sections/zero_chem_pot/data

cd calculations/su3_3d_cutoff_int_cross_sections/zero_chem_pot/data

filename="integratedCrossSections_setA_TMin0p120000_TMax0p300000.ini"
sed -i 's/approximationMethod = COMPLETE_COV/approximationMethod = ZHUANG/g' $filename
./nambuJonaLasinioModel.out use-config-file $filename
sed -i 's/approximationMethod = ZHUANG/approximationMethod = COMPLETE_COV/g' $filename

filename="integratedCrossSections_setA_TMin0p205000_TMax0p220000.ini"
sed -i 's/approximationMethod = COMPLETE_COV/approximationMethod = ZHUANG/g' $filename
./nambuJonaLasinioModel.out use-config-file $filename
sed -i 's/approximationMethod = ZHUANG/approximationMethod = COMPLETE_COV/g' $filename

rm nambuJonaLasinioModel.out

cd ..
