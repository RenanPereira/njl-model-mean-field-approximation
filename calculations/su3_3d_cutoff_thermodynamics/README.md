# README

## Overview

This folder contains the necessary scripts, data files, and configuration files to evaluate the effective masses and thermodynamics of the SU(3) NJL model for different values of temperature and chemical potential using different parameter sets.

Regarding the parameter sets, we consider sets both without (Set A) and with 8-quark interactions at the Lagrangian level (Sets B and C). For more details about the different parameter sets, see [here](../su3_3d_cutoff_phase_diagram/README.md).

### Folder Structure
```
.
├── README.md                # This file
├── build_plots.sh           # Shell script to build plots
├── __init__.py              # Python file necessary to modularize
├── data                     # Contains input .ini files and generated data files
│   ├── *.ini                # Configuration files for the calculations
│   ├── *.dat                # Generated data files
|   ├── main.cpp             # Local `main.cpp` file with hardcoded scenarios for testing
├── execute_calculations.sh  # Shell script to execute the calculations
├── execute_local_main.sh    # Shell script that swaps the `main.cpp` file of the source code by a local one,
|                              compiles and runs the executable
├── plotting                 # Python scripts for generating plots
│   ├── build_plots_*.py     # Specific plot scripts for various scenarios
│   ├── __init__.py          # Python file necessary to modularize
├── plots                    # Directory to store the generated plot images
│   ├── *.png                # Output plot files
```

### Prerequisites
- Ensure the necessary build tools and compilers are installed (e.g., `make` and a C++ compiler).
- Python 3.x must be installed along with the required libraries for plotting (e.g., `matplotlib`, `numpy`).

## How to Use

### Step 1: Execute Calculations
Run the scripts inside the folders (e.g., `fixed_chem_pot_temp`) to compile the C++ code, execute it using the provided configuration files, and generate the data files. For instance, inside the folder `fixed_chem_pot_temp`, you can find:

```bash
./execute_calculations.sh
./execute_local_main.sh
```
These scripts:
1. Navigates to the parent directory and builds the C++ code using `make`.
2. Copies the compiled binary to the appropriate folder.
3. Executes the binary with different `.ini` configuration files (or, in the case of `execute_local_main.sh`, compiles the code using a different `main.cpp` file).
4. Cleans up the binary after execution.

### Step 2: Generate Plots
Run the `build_plots.sh` script to generate the plots based on the generated data files.

```bash
./build_plots.sh
```
This script:
1. Navigates to the previous directory.
2. Executes the Python scripts to build the plots for each scenario.

### Running the local `main.cpp` file:
As described before, there is a script that switches the `main.cpp` file from the root project by the `main.cpp` file in the `\data` folder before compiling the code.
```bash
./execute_local_main.sh
```

## Output
- **Data files**: Stored in the `data` directory, with `.dat` extensions.
- **Plots**: Saved in the `plots` directory, with `.png` extensions.

## Notes
- Ensure that the `plots` and `data` directories are writable by the scripts.


### Example Commands

To execute everything in one go:
```bash
cd fixed_chem_pot_temp
./execute_calculations.sh
./build_plots.sh
```

## Results Zero chemical potential

In this section we present some thermodynamic quantities of the SU3 NJL model resulting from these calculations.

### Quark masses

<p align="center">
  <img src="fixed_chem_pot_temp/plots/quark_eff_masses_vs_temp_CP0_setA.png" width="32%">
  <img src="fixed_chem_pot_temp/plots/quark_eff_masses_vs_temp_CP0_setB.png" width="32%">
  <img src="fixed_chem_pot_temp/plots/quark_eff_masses_vs_temp_CP0_setC.png" width="32%">
</p>

<p align="center">
  <img src="fixed_chem_pot_temp/plots/quark_eff_masses_vs_temp_CP0_setsABC.png" width="32%">
</p>

### Pressure

<p align="center">
  <img src="fixed_chem_pot_temp/plots/pressure_vs_temp_CP0_setA.png" width="32%">
  <img src="fixed_chem_pot_temp/plots/pressure_vs_temp_CP0_setB.png" width="32%">
  <img src="fixed_chem_pot_temp/plots/pressure_vs_temp_CP0_setC.png" width="32%">
</p>

<p align="center">
  <img src="fixed_chem_pot_temp/plots/pressure_vs_temp_CP0_setsABC.png" width="32%">
</p>



### Entropy density

<p align="center">
  <img src="fixed_chem_pot_temp/plots/entropy_vs_temp_CP0_setA.png" width="32%">
  <img src="fixed_chem_pot_temp/plots/entropy_vs_temp_CP0_setB.png" width="32%">
  <img src="fixed_chem_pot_temp/plots/entropy_vs_temp_CP0_setC.png" width="32%">
</p>

<p align="center">
  <img src="fixed_chem_pot_temp/plots/entropy_vs_temp_CP0_setsABC.png" width="32%">
  <img src="fixed_chem_pot_temp/plots/s_over_temp3_vs_temp_CP0_setsABC.png" width="32%">
</p>


In the plots below we also calculate the temperature derivative of the pressure at fixed chemical potential, which yields the entropy density. Hence, one can assess the compatibility of the results provided in the data file.

<p align="center">
  <img src="fixed_chem_pot_temp/plots/entropy_dPdT_vs_temp_CP0_setA.png" width="32%">
  <img src="fixed_chem_pot_temp/plots/entropy_dPdT_vs_temp_CP0_setB.png" width="32%">
  <img src="fixed_chem_pot_temp/plots/entropy_dPdT_vs_temp_CP0_setC.png" width="32%">
</p>

### Energy density

<p align="center">
  <img src="fixed_chem_pot_temp/plots/energy_vs_temp_CP0_setA.png" width="32%">
  <img src="fixed_chem_pot_temp/plots/energy_vs_temp_CP0_setB.png" width="32%">
  <img src="fixed_chem_pot_temp/plots/energy_vs_temp_CP0_setC.png" width="32%">
</p>

<p align="center">
  <img src="fixed_chem_pot_temp/plots/energy_vs_temp_CP0_setsABC.png" width="32%">
</p>

In the plots below we also calculate the energy density via the Euler equation. Hence, one can assess the compatibility of the results provided in the data file.

<p align="center">
  <img src="fixed_chem_pot_temp/plots/energy_euler_eq_vs_temp_CP0_setA.png" width="32%">
  <img src="fixed_chem_pot_temp/plots/energy_euler_eq_vs_temp_CP0_setB.png" width="32%">
  <img src="fixed_chem_pot_temp/plots/energy_euler_eq_vs_temp_CP0_setC.png" width="32%">
</p>

### Pressure and Energy density

<p align="center">
  <img src="fixed_chem_pot_temp/plots/pressure_vs_energy_CP0_setA.png" width="32%">
  <img src="fixed_chem_pot_temp/plots/pressure_vs_energy_CP0_setB.png" width="32%">
  <img src="fixed_chem_pot_temp/plots/pressure_vs_energy_CP0_setC.png" width="32%">
</p>

<p align="center">
  <img src="fixed_chem_pot_temp/plots/pressure_vs_energy_CP0_setsABC.png" width="32%">
</p>
