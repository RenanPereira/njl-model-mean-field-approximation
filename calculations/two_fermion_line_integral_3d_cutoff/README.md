# README

## Overview
This folder contains the necessary scripts, data files, and configuration files to evaluate the two fermion line integral for various scenarios of momentum, temperature, and chemical potential. It also includes tools to generate corresponding plots for visualization.

The majority of these results have been published in the following paper:
[A new approach to the 3-momentum regularization of the in-medium one and two fermion line integrals with applications to cross sections in the Nambu–Jona-Lasinio model](https://arxiv.org/pdf/2310.05749)
Or, check the [DOI](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.109.025206).


### Folder Structure
```
.
├── README.md                # This file
├── build_plots.sh           # Shell script to build plots
├── __init__.py              # Python file necessary to modularize
├── data                     # Contains input .ini files and generated data files
│   ├── *.ini                # Configuration files for the calculations
│   ├── *.dat                # Generated data files
├── execute_calculations.sh  # Shell script to execute the calculations
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

## Output
- **Data files**: Stored in the `data` directory, with `.dat` extensions.
- **Plots**: Saved in the `plots` directory, with `.png` extensions.

## Notes
- Modify the `.ini` files in the `data` directory to configure specific scenarios of momentum, temperature, and chemical potential.
- Ensure that the `plots` and `data` directories are writable by the scripts.

### Example Commands
To execute everything in one go:
```bash
./execute_calculations.sh
./build_plots.sh
```

