#!/bin/bash

# This script runs all the calculations defined inside the calculations folder.
# This can be used to test the code base and understand if the modification of the 
# code or implementation of new features broke something unexpectedly.

cd ../../


# B0 Integral Study
echo "Executing: B0 Integral Study"
(cd calculations/two_fermion_line_integral_3d_cutoff && ./execute_calculations.sh)


# SU3 NJL Vacuum Study with and without 8q interactions
echo "Executing: SU3 NJL Vacuum Study with and without 8q interactions"
(cd calculations/su3_3d_cutoff_vacuum_masses && ./execute_calculations.sh)


# SU3 NJL Phase Diagram Study with and without 8q interactions
echo "Executing: SU3 NJL Phase Diagram Study with and without 8q interactions"
(cd calculations/su3_3d_cutoff_phase_diagram && ./execute_calculations.sh)


# SU3 NJL Cross Section Study with Klevansky parameter set
# This calculation uses parallel computing for better permoance. One has to adapt the calculation
# to the number of threads in the machine running the test.
echo "Executing: SU3 NJL Cross Section Study with Klevansky parameter set"
threads=$(nproc)
echo "This machine has $threads threads."
threads=$((threads - 1))
echo "The calculation will be executed using $threads threads."
sed -i "s/^numberOfThreads = .*/numberOfThreads = $threads/" calculations/su3_3d_cutoff_cross_sections_klevansky/data/crossSections_T0.215000_CP0.000000.ini
sed -i "s/^numberOfThreads = .*/numberOfThreads = $threads/" calculations/su3_3d_cutoff_cross_sections_klevansky/data/crossSections_T0.250000_CP0.000000.ini
sed -i "s/^numberOfThreads = .*/numberOfThreads = $threads/" calculations/su3_3d_cutoff_cross_sections_klevansky/data/crossSections_T0.250000_CP0.100000.ini
sed -i "s/^numberOfThreads = .*/numberOfThreads = $threads/" calculations/su3_3d_cutoff_cross_sections_klevansky/data/crossSections_T0.250000_CP0.200000.ini
sed -i "s/^numberOfThreads = .*/numberOfThreads = $threads/" calculations/su3_3d_cutoff_cross_sections_klevansky/data/crossSections_T0.250000_CP0.300000.ini
sed -i "s/^numberOfThreads = .*/numberOfThreads = $threads/" calculations/su3_3d_cutoff_cross_sections_klevansky/data/crossSections_T0.250000_CP0.400000.ini
sed -i "s/^numberOfThreads = .*/numberOfThreads = $threads/" calculations/su3_3d_cutoff_cross_sections_klevansky/data/crossSections_T0.250000_CP0.500000.ini
(cd calculations/su3_3d_cutoff_cross_sections_klevansky && ./execute_calculations.sh)


# SU3 NJL Integrated Cross Section Study
echo "Executing: SU3 NJL Integrated Cross Section Study"
(cd calculations/su3_3d_cutoff_int_cross_sections/zero_chem_pot && ./execute_calculations_COMPLETE_COV.sh)
(cd calculations/su3_3d_cutoff_int_cross_sections/zero_chem_pot && ./execute_calculations_KLEVANSKY.sh)
(cd calculations/su3_3d_cutoff_int_cross_sections/zero_chem_pot && ./execute_calculations_ZHUANG.sh)


# SU3 NJL Quark Relaxation Time Study
echo "Executing: SU3 NJL Quark Relaxation Time Study"
(cd calculations/su3_3d_cutoff_quark_relaxation_times && ./execute_calculations.sh)


# SU3 NJL Transport Coefficients Study
echo "Executing: SU3 NJL Transport Coefficients Study"
(cd calculations/su3_3d_cutoff_transport_coefficients && ./execute_calculations.sh)


# SU3 NJL Thermodynamics Study with and without 8q interactions
echo "Executing: SU3 NJL Thermodynamics Study with and without 8q interactions"
(cd calculations/su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp && ./execute_calculations.sh)
