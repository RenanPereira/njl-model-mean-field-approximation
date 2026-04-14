#!/bin/bash

echo "Script that computes the quark relaxation time within the NJL model using the data from the integrated cross section study"

cd ..

python3 -m su3_3d_cutoff_quark_relaxation_times.compute_quark_rel_times_cp0
python3 -m su3_3d_cutoff_quark_relaxation_times.compute_quark_rel_times_cpcep
