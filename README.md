This file is a work in progress...

Compile main code:
```
make
# or, using parallel computing,
make -j$(nproc)
```
to clean the artifacts, use:
```
make clean
```

To test the main code by including a `.ini` configuration file:
```
(cd bin/ && ./nambuJonaLasinioModel.out use-config-file tests.ini)
```

To test the `ini_file_parser` module, one can execute the `execute_tests.sh` script present in that folder. For that, use:
```
(cd src/ini_file_parser/ && ./execute_tests.sh)
```

To test the `integration_methods` module, one can execute the `execute_tests.sh` script present in that folder. For that, use:
```
(cd src/integration_methods/ && ./execute_tests.sh)
```

To test the `gsl_wrapper` module, one can execute the `execute_tests.sh` script present in that folder. For that, use:
```
(cd src/gsl_wrapper/ && ./execute_tests.sh)
```

# Calculations

## B0 Integral Study (two_fermion_line_integral_3d_cutoff)
```
(cd calculations/two_fermion_line_integral_3d_cutoff && ./execute_calculations.sh)
(cd calculations/two_fermion_line_integral_3d_cutoff && ./build_plots.sh)
```

## SU3 NJL Vacuum (su3_3d_cutoff_vacuum_masses)

### Vacuum Masses with and without 8q interactions
```
(cd calculations/su3_3d_cutoff_vacuum_masses && ./execute_calculations.sh)
```

## SU3 NJL Phase Diagram Study (su3_3d_cutoff_phase_diagram)
```
(cd calculations/su3_3d_cutoff_phase_diagram && ./execute_calculations.sh)
(cd calculations/su3_3d_cutoff_phase_diagram && ./build_plots.sh)
```

## SU3 NJL Cross Section Study 

### Klevansky parameter set (su3_3d_cutoff_cross_sections_klevansky)
```
(cd calculations/su3_3d_cutoff_cross_sections_klevansky && ./execute_calculations.sh)
(cd calculations/su3_3d_cutoff_cross_sections_klevansky && ./build_plots.sh)
```

## SU3 NJL Integrated Cross Section Study (su3_3d_cutoff_integrated_cross_sections)
```
(cd calculations/su3_3d_cutoff_int_cross_sections/zero_chem_pot && ./execute_calculations_COMPLETE_COV.sh)
(cd calculations/su3_3d_cutoff_int_cross_sections/zero_chem_pot && ./execute_calculations_KLEVANSKY.sh)
(cd calculations/su3_3d_cutoff_int_cross_sections/zero_chem_pot && ./execute_calculations_ZHUANG.sh)
(cd calculations/su3_3d_cutoff_int_cross_sections/zero_chem_pot && ./execute_local_main.sh)
(cd calculations/su3_3d_cutoff_int_cross_sections/zero_chem_pot && ./build_plots.sh)
```

## SU3 NJL Quark Relaxation Time Study (su3_3d_cutoff_quark_relaxation_times)
```
(cd calculations/su3_3d_cutoff_quark_relaxation_times && ./execute_calculations.sh)
(cd calculations/su3_3d_cutoff_quark_relaxation_times && ./build_plots.sh)
```


# Conventions


This code was written by a Physicist without the most in-depth knowledge about the standards and structures of `Clean Code`. The following conventions are being used in the code (at least trying to...): 

C++ code:
- classes  → PascalCase
- namespaces → PascalCase
- methods → camelCase
- member variables → camelCase
- header guard → UPPER_CASE with underscores

Python code:
- 
