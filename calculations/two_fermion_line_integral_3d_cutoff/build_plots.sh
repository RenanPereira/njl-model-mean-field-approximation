#!/bin/bash

echo "Script that builds plots of the two fermion line integral for various scenarios"
echo ""

cd ..

python3 -m two_fermion_line_integral_3d_cutoff.plotting.build_plots_B0_vs_k0
python3 -m two_fermion_line_integral_3d_cutoff.plotting.build_plots_B0_vs_k
python3 -m two_fermion_line_integral_3d_cutoff.plotting.build_plots_B0_vs_k0_diff_T
python3 -m two_fermion_line_integral_3d_cutoff.plotting.build_plots_B0_vs_k_diff_T
python3 -m two_fermion_line_integral_3d_cutoff.plotting.build_plots_B0_vs_k0_diff_mu
python3 -m two_fermion_line_integral_3d_cutoff.plotting.build_plots_B0_vs_k_diff_mu
python3 -m two_fermion_line_integral_3d_cutoff.plotting.build_plots_B0_vs_k0_diff_T_mu
python3 -m two_fermion_line_integral_3d_cutoff.plotting.build_plots_B0_vs_k_diff_T_mu

cd ..
