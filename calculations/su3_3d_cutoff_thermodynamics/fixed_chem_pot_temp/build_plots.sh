#!/bin/bash

echo "Script that builds plots of the NJL model thermodynamics"
echo ""

cd ../../

python3 -m su3_3d_cutoff_thermodynamics.fixed_chem_pot_temp.plotting.generate_effective_masses
python3 -m su3_3d_cutoff_thermodynamics.fixed_chem_pot_temp.plotting.generate_entropy_density
python3 -m su3_3d_cutoff_thermodynamics.fixed_chem_pot_temp.plotting.generate_pressure
python3 -m su3_3d_cutoff_thermodynamics.fixed_chem_pot_temp.plotting.generate_energy_density_plots
python3 -m su3_3d_cutoff_thermodynamics.fixed_chem_pot_temp.plotting.generate_pressure_energy_density
python3 -m su3_3d_cutoff_thermodynamics.fixed_chem_pot_temp.plotting.generate_quark_density_plots
