#!/bin/bash

# This script runs a portion of the calculations defined inside the calculations folder.
# This can be used to test the code base and understand if the modification of the 
# code or implementation of new features broke something unexpectedly.

cd ../.. || exit

# B0 Integral Study
(cd calculations/two_fermion_line_integral_3d_cutoff && ./execute_calculations.sh)

# SU3 NJL Vacuum Study with and without 8q interactions
(cd calculations/su3_3d_cutoff_vacuum_masses && ./execute_calculations.sh)

# SU3 NJL Phase Diagram Study with and without 8q interactions
(cd calculations/su3_3d_cutoff_phase_diagram && ./execute_calculations.sh)

# SU3 NJL Cross Section Study with Klevansky parameter set
(cd calculations/su3_3d_cutoff_cross_sections_klevansky && ./execute_calculations.sh)

# SU3 NJL Thermodynamics Study with and without 8q interactions
(cd calculations/su3_3d_cutoff_thermodynamics/fixed_chem_pot_temp && ./execute_calculations.sh)

# SU3 NJL Integrated Cross Section Study
(cd calculations/su3_3d_cutoff_int_cross_sections/zero_chem_pot && ./execute_calculations_COMPLETE_COV.sh)

# SU3 NJL Quark Relaxation Time Study
(cd calculations/su3_3d_cutoff_quark_relaxation_times && ./execute_calculations.sh)

# SU3 NJL Transport Coefficients Study
(cd calculations/su3_3d_cutoff_transport_coefficients && ./execute_calculations.sh)
