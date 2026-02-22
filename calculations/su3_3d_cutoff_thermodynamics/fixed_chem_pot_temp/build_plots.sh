#!/bin/bash

echo "Script that builds plots of the NJL model thermodynamics"
echo ""

cd ../../

python3 -m su3_3d_cutoff_thermodynamics.fixed_chem_pot_temp.plotting.build_plots_thermodynamics
