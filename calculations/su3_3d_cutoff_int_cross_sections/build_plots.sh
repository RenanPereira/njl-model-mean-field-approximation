#!/bin/bash

echo "Script that builds plots of the NJL model integrated cross sections at the CEP chemical potential and using different methods to evaluate the integral"
echo ""

cd ../

python3 -m su3_3d_cutoff_int_cross_sections.plotting.generate_plots_int_cross_sections_setA_CP0p000000
python3 -m su3_3d_cutoff_int_cross_sections.plotting.generate_plots_int_cross_sections_setA_CP0p318436
