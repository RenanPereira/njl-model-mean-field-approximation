#!/bin/bash

echo "Script that computes transport coefficients within the NJL model using the data from the quark relaxation times study"

cd ..

python3 -m su3_3d_cutoff_transport_coefficients.compute_shear_viscosity
