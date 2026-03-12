#!/bin/bash

echo "Script that calculates the Cross Sections for the NJL model using different values of temperature and chemical potential"
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

make -j"$threads"

cp bin/nambuJonaLasinioModel.out calculations/su3_3d_cutoff_cross_sections_klevansky/data

cd calculations/su3_3d_cutoff_cross_sections_klevansky/data

sed -i "s/^numberOfThreads = .*/numberOfThreads = $threads/" crossSections_T0p215000_CP0p000000.ini
sed -i "s/^numberOfThreads = .*/numberOfThreads = $threads/" crossSections_T0p250000_CP0p000000.ini
sed -i "s/^numberOfThreads = .*/numberOfThreads = $threads/" crossSections_T0p250000_CP0p100000.ini
sed -i "s/^numberOfThreads = .*/numberOfThreads = $threads/" crossSections_T0p250000_CP0p200000.ini
sed -i "s/^numberOfThreads = .*/numberOfThreads = $threads/" crossSections_T0p250000_CP0p300000.ini
sed -i "s/^numberOfThreads = .*/numberOfThreads = $threads/" crossSections_T0p250000_CP0p400000.ini
sed -i "s/^numberOfThreads = .*/numberOfThreads = $threads/" crossSections_T0p250000_CP0p500000.ini

./nambuJonaLasinioModel.out use-config-file crossSections_T0p215000_CP0p000000.ini
./nambuJonaLasinioModel.out use-config-file crossSections_T0p250000_CP0p000000.ini
./nambuJonaLasinioModel.out use-config-file crossSections_T0p250000_CP0p100000.ini
./nambuJonaLasinioModel.out use-config-file crossSections_T0p250000_CP0p200000.ini
./nambuJonaLasinioModel.out use-config-file crossSections_T0p250000_CP0p300000.ini
./nambuJonaLasinioModel.out use-config-file crossSections_T0p250000_CP0p400000.ini
./nambuJonaLasinioModel.out use-config-file crossSections_T0p250000_CP0p500000.ini

rm nambuJonaLasinioModel.out

sed -i "s/^numberOfThreads = .*/numberOfThreads = 15/" crossSections_T0p215000_CP0p000000.ini
sed -i "s/^numberOfThreads = .*/numberOfThreads = 15/" crossSections_T0p250000_CP0p000000.ini
sed -i "s/^numberOfThreads = .*/numberOfThreads = 15/" crossSections_T0p250000_CP0p100000.ini
sed -i "s/^numberOfThreads = .*/numberOfThreads = 15/" crossSections_T0p250000_CP0p200000.ini
sed -i "s/^numberOfThreads = .*/numberOfThreads = 15/" crossSections_T0p250000_CP0p300000.ini
sed -i "s/^numberOfThreads = .*/numberOfThreads = 15/" crossSections_T0p250000_CP0p400000.ini
sed -i "s/^numberOfThreads = .*/numberOfThreads = 15/" crossSections_T0p250000_CP0p500000.ini

cd ..
