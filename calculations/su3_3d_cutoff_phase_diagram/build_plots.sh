#!/bin/bash

echo "Script that builds plots of the NJL model phase diagram for different parameter sets"
echo ""

cd ..

python3 -m su3_3d_cutoff_phase_diagram.plots_scripts.build_plots_phase_diagram_chemPot_vs_temp
python3 -m su3_3d_cutoff_phase_diagram.plots_scripts.build_plots_phase_diagram_rhoB_vs_temp
