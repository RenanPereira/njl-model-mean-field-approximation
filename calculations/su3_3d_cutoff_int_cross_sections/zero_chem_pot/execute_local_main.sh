# This script must be executed from the root folder using: 
# cd calculations/su3_3d_cutoff_integrated_cross_sections && ./execute_local_main.sh && cd ../../

#!/bin/bash

echo "Script that switches the main.cpp file from the root project by the main.cpp file in this folder before compiling the code."
echo "In the local main.cpp file, Integrated Cross Sections for the NJL model are calculated for different parameter sets \
as functions of temperature and chemical potential."
echo ""

# Go to project root
cd .. && cd .. && cd ..

# Clean previous build 
make clean

# Backup original main.cpp file for later restoration
cp src/main.cpp src/main_original.cpp && rm src/main.cpp

# Copy local main.cpp to root folder of the project for build step
cp calculations/su3_3d_cutoff_int_cross_sections/zero_chem_pot/data/main.cpp src/main.cpp

# Compile code with the local main.cpp file
make -j$(nproc)

# Restore to original main.cpp file
rm src/main.cpp && cp src/main_original.cpp src/main.cpp && rm src/main_original.cpp

# Copy executable file to data folder where it will be executed
cp bin/nambuJonaLasinioModel.out calculations/su3_3d_cutoff_int_cross_sections/zero_chem_pot/data

# Run executable
cd calculations/su3_3d_cutoff_int_cross_sections/zero_chem_pot/data
./nambuJonaLasinioModel.out
rm nambuJonaLasinioModel.out

# Go to project root
cd ..
