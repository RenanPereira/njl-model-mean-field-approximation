#!/bin/bash

echo "Script that calculates the Integrated Cross Sections for the NJL model."
echo "This calculation is for testing purposes only: a small sample of points are calculated."
echo ""

cd ../.. || exit

# Set threads: use argument if provided, otherwise use all cores
echo "This calculation can be parallelized."
machine_threads=$(nproc)
echo "This machine has a total of $machine_threads threads."

if [ -n "$1" ]; then
    if [[ "$1" =~ ^[0-9]+$ ]] && [ "$1" -le "$machine_threads" ]; then
        threads=$1
    else
        echo "Error: Please provide a valid number of threads ≤ $machine_threads."
        exit 1
    fi
else
    threads=$((machine_threads - 1))
fi
echo "The calculation will be executed using $threads threads."
echo ""

# Clean previous build and re-build
make clean
make -j$(nproc)

cp bin/nambuJonaLasinioModel.out calculations/su3_3d_cutoff_int_cross_sections_test/data

filename="integratedCrossSections_setA_TMin0p290000_TMax0p300000.ini"

sed -i "s/^numberOfThreads = .*/numberOfThreads = $threads/" $filename

cd calculations/su3_3d_cutoff_int_cross_sections_test/data

# setA
./nambuJonaLasinioModel.out use-config-file $filename

rm nambuJonaLasinioModel.out

sed -i "s/^numberOfThreads = .*/numberOfThreads = 15/" $filename

cd ..
