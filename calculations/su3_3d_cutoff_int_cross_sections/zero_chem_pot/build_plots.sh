#!/bin/bash

echo "Script that builds plots of the NJL model integrated cross sections for zero chemical potential and using different methods to evaluate the integral"
echo ""

cd ../../

python3 -m su3_3d_cutoff_int_cross_sections.zero_chem_pot.plotting.build_plots_integrated_cross_sections
