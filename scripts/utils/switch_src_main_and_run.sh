#!/bin/bash

# This script automates the process of:
# - switching the main file in the src folder by a local one present at the path provided by the user;
# - building the entire C++ code using it the provided file;
# - copies the executable to the same path as the one provided by the user;
# - execute the code.
#
# It takes three arguments: the path to the main file which will be used to substitute the 
# one in the src directory, the name of the file and the number of folder to cd up from 
# the directory where the script is being executed.
# 
# The script creates a backup of the original file main file which is restored after the code is builded
# with the main file present at the path provided by the user.

# Validate the input arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <main_folder_path> <local_main_filename> <number_of_up_folders_to_makefile>"
    exit 1
fi

main_folder_path="$1"
local_main_filename="$2"
number_of_up_folders_to_makefile="$3"

echo "Running the script that switches the main.cpp file from the src folder by the main.cpp file in the folder provided."
echo ""

# Go to project root where Makefile is found
for ((i=0; i<number_of_up_folders_to_makefile; i++)); do
  cd ..
done

# Clean previous build 
make clean

# Backup original main.cpp file for later restoration
cp src/main.cpp src/main.bak && rm src/main.cpp

# Copy local main.cpp to root folder of the project for build step
cp "${main_folder_path}"/"${local_main_filename}" src/main.cpp

# Compile code with the local main.cpp file
make -j$(nproc)

# Restore to original main.cpp file
rm src/main.cpp && cp src/main.bak src/main.cpp && rm src/main.bak

# Copy executable file to the folder where it will be executed
cp bin/nambuJonaLasinioModel.out "${main_folder_path}"

# Run executable
cd "$main_folder_path"
./nambuJonaLasinioModel.out
rm nambuJonaLasinioModel.out

# Go to project root
cd ..
