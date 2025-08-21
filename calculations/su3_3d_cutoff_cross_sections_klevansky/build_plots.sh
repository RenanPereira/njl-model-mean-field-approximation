#!/bin/bash

echo "Script that generates plots of cross sections in the NJL model"
echo ""

cd ..

python3 -m su3_3d_cutoff_cross_sections_klevansky.plots_scripts.build_plots_cross_sections_zero_chem_pot
python3 -m su3_3d_cutoff_cross_sections_klevansky.plots_scripts.build_plots_cross_sections_finite_chem_pot
