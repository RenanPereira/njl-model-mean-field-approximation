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

## Results

In this section we present the results of numerical results for the two fermion line integral, $B_0$, by choosing a particular value for the cutoff, $\Lambda$, and constant values for the fermion masses, $M_i$ and $M_j$. 


### Different pole shifts

#### External momentum shift vs mass shift
![B0 plot](plots/B0_vs_k0_T0.0Cpi0.0Cpj0.0L1.0Mi0.4Mj0.4k0.0_k0Shift.png)
![B0 plot](plots/B0_vs_k0_T0.0Cpi0.0Cpj0.0L1.0Mi0.4Mj0.4k0.0.png)


#### $B_0$ vs $|\bold{k}|$ different $k_0$

![B0 plot](plots/B0_vs_k_T0.0Cpi0.0Cpj0.0L1.0Mi0.4Mj0.4_diff_k0.png)


### $B_0$ vs $|\bold{k}|$ different $k_0$ and temperature

![B0 plot](plots/B0_vs_k_Cpi0.0Cpj0.0L1.0Mi0.4Mj0.4k00.0_diff_T.png)
![B0 plot](plots/B0_vs_k_Cpi0.0Cpj0.0L1.0Mi0.4Mj0.4k00.5_diff_T.png)
![B0 plot](plots/B0_vs_k_Cpi0.0Cpj0.0L1.0Mi0.4Mj0.4k01.0_diff_T.png)
![B0 plot](plots/B0_vs_k_Cpi0.0Cpj0.0L1.0Mi0.4Mj0.4k01.5_diff_T.png)
![B0 plot](plots/B0_vs_k_Cpi0.0Cpj0.0L1.0Mi0.4Mj0.4k02.0_diff_T.png)
![B0 plot](plots/B0_vs_k_Cpi0.0Cpj0.0L1.0Mi0.4Mj0.4k02.5_diff_T.png)


### $B_0$ vs $|\bold{k}|$ different $k_0$ and chemical potential

![B0 plot](plots/B0_vs_k_T0.0L1.0Mi0.4Mj0.4k00.0_diff_mu.png)
![B0 plot](plots/B0_vs_k_T0.0L1.0Mi0.4Mj0.4k00.5_diff_mu.png)
![B0 plot](plots/B0_vs_k_T0.0L1.0Mi0.4Mj0.4k01.0_diff_mu.png)
![B0 plot](plots/B0_vs_k_T0.0L1.0Mi0.4Mj0.4k01.5_diff_mu.png)
![B0 plot](plots/B0_vs_k_T0.0L1.0Mi0.4Mj0.4k02.0_diff_mu.png)
![B0 plot](plots/B0_vs_k_T0.0L1.0Mi0.4Mj0.4k02.5_diff_mu.png)


### $B_0$ vs $|\bold{k}|$ different $k_0$, temperature and chemical potential

![B0 plot](plots/B0_vs_k_L1.0Mi0.4Mj0.4k00.0_diff_T_mu.png)
![B0 plot](plots/B0_vs_k_L1.0Mi0.4Mj0.4k00.5_diff_T_mu.png)
![B0 plot](plots/B0_vs_k_L1.0Mi0.4Mj0.4k01.0_diff_T_mu.png)
![B0 plot](plots/B0_vs_k_L1.0Mi0.4Mj0.4k01.5_diff_T_mu.png)
![B0 plot](plots/B0_vs_k_L1.0Mi0.4Mj0.4k02.0_diff_T_mu.png)
![B0 plot](plots/B0_vs_k_L1.0Mi0.4Mj0.4k02.5_diff_T_mu.png)


#### $B_0$ vs $k_0$ different $|\bold{k}|$

![B0 plot](plots/B0_vs_k0_T0.0Cpi0.0Cpj0.0L1.0Mi0.4Mj0.4_diff_k.png)

### $B_0$ vs $k_0$ different $|\bold{k}|$ and temperature

![B0 plot](plots/B0_vs_k0_Cpi0.0Cpj0.0L1.0Mi0.4Mj0.4k0.0_diff_T.png)
![B0 plot](plots/B0_vs_k0_Cpi0.0Cpj0.0L1.0Mi0.4Mj0.4k0.5_diff_T.png)
![B0 plot](plots/B0_vs_k0_Cpi0.0Cpj0.0L1.0Mi0.4Mj0.4k1.0_diff_T.png)
![B0 plot](plots/B0_vs_k0_Cpi0.0Cpj0.0L1.0Mi0.4Mj0.4k1.5_diff_T.png)

### $B_0$ vs $k_0$ different $|\bold{k}|$ and chemical potential

![B0 plot](plots/B0_vs_k0_T0.0L1.0Mi0.4Mj0.4k0.0_diff_mu.png)
![B0 plot](plots/B0_vs_k0_T0.0L1.0Mi0.4Mj0.4k0.5_diff_mu.png)
![B0 plot](plots/B0_vs_k0_T0.0L1.0Mi0.4Mj0.4k1.0_diff_mu.png)
![B0 plot](plots/B0_vs_k0_T0.0L1.0Mi0.4Mj0.4k1.5_diff_mu.png)

### $B_0$ vs $k_0$ different $|\bold{k}|$, temperature and chemical potential

![B0 plot](plots/B0_vs_k0_L1.0Mi0.4Mj0.4k0.0_diff_T_mu.png)
![B0 plot](plots/B0_vs_k0_L1.0Mi0.4Mj0.4k0.5_diff_T_mu.png)
![B0 plot](plots/B0_vs_k0_L1.0Mi0.4Mj0.4k1.0_diff_T_mu.png)
![B0 plot](plots/B0_vs_k0_L1.0Mi0.4Mj0.4k1.5_diff_T_mu.png)
