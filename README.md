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
cd bin/ && ./nambuJonaLasinioModel.out use-config-file tests.ini && cd ..
```

To test the `ini_file_parser` module, one can execute the `execute_tests.sh` script present in that folder. For that, use:
```
cd src/ini_file_parser/ && ./execute_tests.sh && cd ../../
```

To test the `integration_methods` module, one can execute the `execute_tests.sh` script present in that folder. For that, use:
```
cd src/integration_methods/ && ./execute_tests.sh && cd ../../
```

To test the `gsl_wrapper` module, one can execute the `execute_tests.sh` script present in that folder. For that, use:
```
cd src/gsl_wrapper/ && ./execute_tests.sh && cd ../../
```
