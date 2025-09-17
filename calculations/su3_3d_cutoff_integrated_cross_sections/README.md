# README

## Overview

This folder contains the necessary scripts, data files, and configuration files to evaluate the integrated cross sections of the SU(3) NJL model for different values of temperature and chemical potential. Additionally, different NJL parameter sets are used, as well as, different "approximations" to evaluate the integral over the differential cross sections.

Regarding the parameter sets, we consider sets both without (Set A) and with 8-quark interactions at the Lagrangian level (Sets B and C). The NJL parameter set A, is the usual Klevansky parameter set. The NJL parameters given in Sets B and C contain 8-quark interactions. In these sets the coupling $g_1$ was fixed manually and the remaining six free parameters were found by requiring the model to reproduce the masses of the pion ($M_{\pi^\pm} = 0.140 \,\mathrm{GeV}$), the kaon ($M_{K^\pm} = 0.494 \,\mathrm{GeV}$), the eta prime ($M_\eta'= 0.958 \,\mathrm{GeV}$) and $a_0^\pm$ ( $M_{a_0^\pm}= 0.960 \,\mathrm{GeV}$) mesons, the leptonic decays of the pion ($f_{\pi^\pm} = 0.0924 \,\mathrm{GeV}$) and kaon ($f_{K^\pm}=0.094 \,\mathrm{GeV}$). For further details, see [Renan Camara Pereira PhD thesis](https://estudogeral.uc.pt/handle/10316/95294).

Regarding the different "approximations" to evaluate the integral over the differential cross sections, 4 values are allowed in the code:
- `COMPLETE_OG`
- `COMPLETE_COV` 
- `KLEVANSKY`
- `ZHUANG`

The first two methods, `COMPLETE_OG` and `COMPLETE_COV`, are a complete evaluation of the integral, i.e., no approximations is made. In the first method however, the 3-dimensional integral is made in its completeness, making the numerical calculation slower, more complex and prone to larger numerical errors due to the precisions one has to take in order to evaluate the integral. The second method, `COMPLETE_COV`, uses a change of variables making nly necessary to evaluate a 2-dimensional integral. Furthermore, the evaluation of the first integral does not involve the evaluation of the differential cross section, making it much more optimized when compared with the `COMPLETE_OG` approach. Of course, both methods must return the same answer, aside from numerical errors. 

The third method, `KLEVANSKY`, is the one employed by Klevansky et. al. in the paper [Elastic scattering and transport coefficients for a quark plasma in SUf(3) at finite temperatures](https://www.sciencedirect.com/science/article/abs/pii/0375947496002473).

The fourth method, `ZHUANG`, is the one suggested by Zhuang et. al. in the paper [Transport properties of a quark plasma and critical scattering at the chiral phase transition](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.51.3728). 

THE DISCUSSION OVER THE "APPROXIMATION" METHODS MUST BE IMPROVED.


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
├── plots_scripts            # Python scripts for generating plots
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
Run the `execute_calculations.sh` script to compile the C++ code, execute it using the provided configuration files, and generate the data files.

```bash
./execute_calculations.sh
```
This script:
1. Navigates to the parent directory and builds the C++ code using `make`.
2. Copies the compiled binary to the appropriate folder.
3. Executes the binary with different `.ini` configuration files.
4. Cleans up the binary after execution.

### Step 2: Generate Plots
Run the `build_plots.sh` script to generate the plots based on the generated data files.

```bash
./build_plots.sh
```
This script:
1. Navigates to the previous directory.
2. Executes the Python scripts to build the plots for each scenario.

### Optional Step: Running the local `main.cpp` file:
Script that switches the `main.cpp` file from the root project by the `main.cpp` file in the `\data` folder before compiling the code.
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
./execute_calculations.sh
./build_plots.sh
```

## Results

In this section we present the results of the integrated cross sections for different processes ($\sigma$ in mb), as a function of temperature and chemical potential.
