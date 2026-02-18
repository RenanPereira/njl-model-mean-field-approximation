#!/bin/bash

echo "Script that builds plots of transport coefficients using different methods to evaluate the integrated cross section and for different parameter sets and physical scenarios"
echo ""

cd ..

python3 -m su3_3d_cutoff_transport_coefficients.plotting.build_plots_shear_viscosity
