#!/bin/bash

# This script builds all the plots defined inside the calculations folder.
# This can be used to test the code base and understand if the modification of the 
# code or implementation of new features broke something unexpectedly.

cd ../../


# B0 Integral Study
echo "Building plots: B0 Integral Study"
(cd calculations/two_fermion_line_integral_3d_cutoff && ./build_plots.sh)


# SU3 NJL Phase Diagram Study with and without 8q interactions
echo "Building plots: SU3 NJL Phase Diagram Study with and without 8q interactions"
(cd calculations/su3_3d_cutoff_phase_diagram && ./build_plots.sh)


# SU3 NJL Cross Section Study with Klevansky parameter set
# This calculation uses parallel computing for better permoance. One has to adapt the calculation
# to the number of threads in the machine running the test.
echo "Building plots: SU3 NJL Cross Section Study with Klevansky parameter set"
(cd calculations/su3_3d_cutoff_cross_sections_klevansky && ./build_plots.sh)


# SU3 NJL Integrated Cross Section Study
echo "Building plots: SU3 NJL Integrated Cross Section Study"
(cd calculations/su3_3d_cutoff_int_cross_sections/zero_chem_pot && ./build_plots.sh)


# SU3 NJL Quark Relaxation Time Study
echo "Building plots: SU3 NJL Quark Relaxation Time Study"
(cd calculations/su3_3d_cutoff_quark_relaxation_times && ./build_plots.sh)


# SU3 NJL Transport Coefficients Study
echo "Building plots: SU3 NJL Transport Coefficients Study"
(cd calculations/su3_3d_cutoff_transport_coefficients && ./build_plots.sh)


# SU3 NJL Thermodynamics Study with and without 8q interactions
echo "Building plots: SU3 NJL Thermodynamics Study with and without 8q interactions"
(cd calculations/su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp && ./build_plots.sh)
