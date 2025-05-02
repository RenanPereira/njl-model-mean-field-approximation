#!/bin/bash

echo "Script that builds plots of the NJL model phase diagram for different parameter sets"
echo ""

cd plots_scripts

python3 build_plots_phase_diagram_chemPot_vs_temp.py
python3 build_plots_phase_diagram_rhoB_vs_temp.py

cd ..
