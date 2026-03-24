#!/bin/bash

echo "Script that builds plots of transport coefficients using different methods to evaluate the integrated cross section and for different parameter sets and physical scenarios"
echo ""

cd ..

python3 -m su3_3d_cutoff_transport_coefficients.plotting.generate_eta_plots
python3 -m su3_3d_cutoff_transport_coefficients.plotting.generate_sigmae_plots
python3 -m su3_3d_cutoff_transport_coefficients.plotting.generate_eta_sigmae_ratios_plots
