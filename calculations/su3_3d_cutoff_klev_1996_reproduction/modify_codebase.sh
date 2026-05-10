#!/bin/bash


echo "Script that modify the codebase to use the 2-line fermion line integral (B0) using \
the Klevansky recipe instead of the usual one defined in two_fermion_line_integral_3d_cutoff.cpp"
echo ""


# Go to project root
cd ../../


# Copy the implementation of the 2-line fermion line integral (B0) using the Klevansky recipe to the proper folder

cp calculations/su3_3d_cutoff_klev_1996_reproduction/two_fermion_line_integral_3d_cutoff_klev_recipe.cpp \
   src/njl_model/n_fermion_line_integrals/two_fermion_line_integral_3d_cutoff_klev_recipe.cpp

cp calculations/su3_3d_cutoff_klev_1996_reproduction/two_fermion_line_integral_3d_cutoff_klev_recipe.h \
   src/njl_model/n_fermion_line_integrals/two_fermion_line_integral_3d_cutoff_klev_recipe.h


# Change the SU3NJL3DCutoffMesonPropagators.cpp file to use the Klevansky recipe

filepath="src/njl_model/su3_3d_cutoff/SU3NJL3DCutoffMesonPropagators.cpp"

# Backup the original file
backup_filepath="${filepath}.bak"
cp "$filepath" "$backup_filepath"

target="two_fermion_line_integral_3d_cutoff.h"
replacement="two_fermion_line_integral_3d_cutoff_klev_recipe.h"
sed -i "s/${target}/${replacement}/g" "$filepath"

target="klevA1 = klevanskyAIntegral3DCutoff(reguScheme, cutoff, T, effCP1, M1, k, integralPrecision)"
replacement="klevA1 = klevanskyAIntegral3DCutoff(reguScheme, cutoff, T, effCP1, M1, 0.0, integralPrecision)"
sed -i "s/${target}/${replacement}/g" "$filepath"

target="klevA2 = klevanskyAIntegral3DCutoff(reguScheme, cutoff, T, effCP2, M2, k, integralPrecision)"
replacement="klevA2 = klevanskyAIntegral3DCutoff(reguScheme, cutoff, T, effCP2, M2, 0.0, integralPrecision)"
sed -i "s/${target}/${replacement}/g" "$filepath"

target="klevB0 = klevanskyB0Integral3DCutoff(reguScheme, T, effCP1, effCP2, cutoff, M1, M2, k0, k, integralPrecision)"
replacement="klevB0 = klevanskyB0Integral3DCutoffKlevanskyRecipe(reguScheme, T, effCP1, effCP2, cutoff, M1, M2, k0, k, integralPrecision)"
sed -i "s/${target}/${replacement}/g" "$filepath"
