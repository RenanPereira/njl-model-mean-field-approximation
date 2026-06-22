# README

## Overview

This folder contains the necessary scripts, data files, and configuration files to evaluate the cross sections of the SU(3) NJL model. This calculation is for testing purposes only: a small sample of points are calculated.


## How to run calculations in this folder

The calculations performed in this folder follow the general structure found in this project. The scripts execute the following steps:

1. Navigates to the parent directory and builds the C++ code using `make`.
2. Copies the compiled binary to the appropriate folder.
3. Executes the binary with different `.ini` configuration files.
4. Cleans up the binary after execution.

The calculations can be executed by running:
```bash
(./execute_calculations.sh)
```
