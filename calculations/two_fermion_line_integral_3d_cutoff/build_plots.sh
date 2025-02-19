#!/bin/bash

echo "Script that builds plots of the two fermion line integral for various scenarios"
echo ""

cd plot_scripts

python3 build_plots_B0_vs_k0.py
python3 build_plots_B0_vs_k.py
python3 build_plots_B0_vs_k0_diff_T.py
python3 build_plots_B0_vs_k_diff_T.py
python3 build_plots_B0_vs_k0_diff_mu.py
python3 build_plots_B0_vs_k_diff_mu.py
python3 build_plots_B0_vs_k0_diff_T_mu.py
python3 build_plots_B0_vs_k_diff_T_mu.py

cd ..
