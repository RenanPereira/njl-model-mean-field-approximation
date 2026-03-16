# Tests

Tests for the `IniFileParser` class are provided in the `main.cpp` file. The tests provided cover all major functionalities:
- `testGetValue()` tests string value retrieval, including non-existent keys.
- `testGetInt()` covers integer retrieval, including cases with missing or malformed data.
- `testGetDouble()` handles double precision values, ensuring correct conversion and handling of defaults.
- `testGetSectionData()` ensures correct section filtering based on prefixes.

In order to execute the tests, simply compile the code in this folder using the Makefile:
```bash
make && ./test.out
make clean
```
In the commands above,
- `make`: Compiles the test suite and the IniFileParser class.
- `./test.out`: Runs the compiled tests.
- `make clean`: Removes the compiled binaries and any temporary files after running the tests.
