#!/bin/bash

echo "Script that builds plots of the quark relaxation time using different methods to evaluate the integrated cross section and for different parameter sets and physical scenarios"
echo ""

cd ..

python3 -m su3_3d_cutoff_quark_relaxation_times.plotting.build_plots_quark_relaxation_times
