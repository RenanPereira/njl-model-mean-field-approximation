#!/bin/bash

# This script automates the process of switching the main file in the src folder by the
# one present at the path provided by the user, builds the entire C++ code using it,
# copies the executable to the same path as the one provided by the user.  
# It takes two arguments: the path to the main file which will be used to substitute the 
# one in the src directory, and the number of folder to cd up from the directory where the 
# script is being executed.
# 
# The script creates a backup of the original file main file which is restored after the code is builded
# with the main file present at the path provided by the user.

# Validate the input arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <main_folder_path> <number_of_up_folders_to_makefile>"
    exit 1
fi

main_folder_path="$1"
number_of_up_folders_to_makefile="$2"

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
cp "${main_folder_path}"/main.cpp src/main.cpp

# Compile code with the local main.cpp file
make -j$(nproc)

# Restore to original main.cpp file
rm src/main.cpp && cp src/main.bak src/main.cpp && rm src/main.bak

# Copy executable file to the folder where it will be executed
cp bin/nambuJonaLasinioModel.out "${main_folder_path}"

# Run executable
cd $main_folder_path
./nambuJonaLasinioModel.out
rm nambuJonaLasinioModel.out

# Go to project root
cd ..
