#!/bin/bash

echo "Script that builds plots of the NJL model thermodynamics"
echo ""

cd ../../

python3 -m su3_3d_cutoff_thermodynamics.fixed_chem_pot_temp.plotting.effective_masses
python3 -m su3_3d_cutoff_thermodynamics.fixed_chem_pot_temp.plotting.entropy_density
python3 -m su3_3d_cutoff_thermodynamics.fixed_chem_pot_temp.plotting.pressure
python3 -m su3_3d_cutoff_thermodynamics.fixed_chem_pot_temp.plotting.energy_density
python3 -m su3_3d_cutoff_thermodynamics.fixed_chem_pot_temp.plotting.pressure_energy_density
python3 -m su3_3d_cutoff_thermodynamics.fixed_chem_pot_temp.plotting.generate_quark_density_plots
