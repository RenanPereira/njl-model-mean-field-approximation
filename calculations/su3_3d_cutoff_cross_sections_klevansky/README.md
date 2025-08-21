# README

## Overview

This folder contains the necessary scripts, data files, and configuration files to evaluate the cross sections of the SU(3) NJL model for different values of temperature and chemical potential. The NJL parameter set used is the usual Klevansky parameter set. For further details, see [New approach to the 3-momentum regularization of the in-medium one- and two-fermion line integrals with applications to cross sections in the Nambu–Jona-Lasinio model](https://arxiv.org/pdf/2310.05749).

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
- Ensure that the `plots` and `data` directories are writable by the scripts.

### Example Commands
To execute everything in one go:
```bash
./execute_calculations.sh
./build_plots.sh
```

## Results

In this section we present the results of the cross sections for different processes ($\sigma$ in mb), as a function of the center-of-mass energy ($\sqrt{s}$ in GeV), taking into account different values of temperature and chemical potential.


### Zero chemical potential

#### uu $\rightarrow$ uu , ud $\rightarrow$ ud (T[GeV]=0.215)
![Cross section plot](plots/cross_section_uuuu_udud_T0215_CP0.png)

#### uu $\rightarrow$ uu , ud $\rightarrow$ ud (T[GeV]=0.250)
![Cross section plot](plots/cross_section_uuuu_udud_T0250_CP0.png)

#### us $\rightarrow$ us , ss $\rightarrow$ ss (T[GeV]=0.215)
![Cross section plot](plots/cross_section_usus_ssss_T0215_CP0.png)

#### us $\rightarrow$ us , ss $\rightarrow$ ss (T[GeV]=0.250)
![Cross section plot](plots/cross_section_usus_ssss_T0250_CP0.png)

#### udbar $\rightarrow$ udbar , uubar $\rightarrow$ uubar, uubar $\rightarrow$ ddbar, uubar $\rightarrow$ ssbar (T[GeV]=0.215)")
![Cross section plot](plots/cross_section_udbarudbar_uubaruubar_uubarddbar_uubarssbar_T0215_CP0.png)

#### udbar $\rightarrow$ udbar , uubar $\rightarrow$ uubar, uubar $\rightarrow$ ddbar, uubar $\rightarrow$ ssbar (T[GeV]=0.250)"
![Cross section plot](plots/cross_section_udbarudbar_uubaruubar_uubarddbar_uubarssbar_T0250_CP0.png)

#### usbar $\rightarrow$ usbar , ssbar $\rightarrow$ uubar, ssbar $\rightarrow$ ssbar (T[GeV]=0.215)
![Cross section plot](plots/cross_section_usbarusbar_ssbaruubar_ssbarssbar_T0215_CP0.png)

#### usbar $\rightarrow$ usbar , ssbar $\rightarrow$ uubar, ssbar $\rightarrow$ ssbar (T[GeV]=0.250)
![Cross section plot](plots/cross_section_usbarusbar_ssbaruubar_ssbarssbar_T0250_CP0.png)


### Finite chemical potential

#### uu $\rightarrow$ uu (T[GeV]=0.250, finite chemical potential)
![Cross section plot](plots/cross_section_uuuu_T0250_CP.png)

#### ubarubar $\rightarrow$ ubarubar (T[GeV]=0.250, finite chemical potential)
![Cross section plot](plots/cross_section_ubarubarubarubar_T0250_CP.png)

#### ud $\rightarrow$ ud (T[GeV]=0.250, finite chemical potential)
![Cross section plot](plots/cross_section_udud_T0250_CP.png)

#### ubardbar $\rightarrow$ ubardbar (T[GeV]=0.250, finite chemical potential)
![Cross section plot](plots/cross_section_ubardbarubardbar_T0250_CP.png)

#### us $\rightarrow$ us (T[GeV]=0.250, finite chemical potential)
![Cross section plot](plots/cross_section_usus_T0250_CP.png)

#### ubarsbar $\rightarrow$ ubarsbar (T[GeV]=0.250, finite chemical potential)
![Cross section plot](plots/cross_section_ubarsbarubarsbar_T0250_CP.png)

#### ss $\rightarrow$ ss (T[GeV]=0.250, finite chemical potential)
![Cross section plot](plots/cross_section_ssss_T0250_CP.png)

#### sbarsbar $\rightarrow$ sbarsbar (T[GeV]=0.250, finite chemical potential)
![Cross section plot](plots/cross_section_sbarsbarsbar_T0250_CP.png)

#### uubar $\rightarrow$ uubar (T[GeV]=0.250, finite chemical potential)
![Cross section plot](plots/cross_section_uubaruubar_T0250_CP.png)

#### uudbar $\rightarrow$ udbar (T[GeV]=0.250, finite chemical potential)
![Cross section plot](plots/cross_section_udbarudbar_T0250_CP.png)

#### uubar $\rightarrow$ ddbar (T[GeV]=0.250, finite chemical potential)
![Cross section plot](plots/cross_section_uubarddbar_T0250_CP.png)

#### uubar $\rightarrow$ ssbar (T[GeV]=0.250, finite chemical potential)
![Cross section plot](plots/cross_section_uubarssbar_T0250_CP.png)

#### ssbar $\rightarrow$ uubar (T[GeV]=0.250, finite chemical potential)
![Cross section plot](plots/cross_section_ssbaruubar_T0250_CP.png)

#### ssbar $\rightarrow$ ssbar (T[GeV]=0.250, finite chemical potential)
![Cross section plot](plots/cross_section_ssbarssbar_T0250_CP.png)

#### usbar $\rightarrow$ usbar (T[GeV]=0.250, finite chemical potential)
![Cross section plot](plots/cross_section_usbarusbar_T0250_CP.png)

#### subar $\rightarrow$ subar (T[GeV]=0.250, finite chemical potential)
![Cross section plot](plots/cross_section_subarsubar_T0250_CP.png)
