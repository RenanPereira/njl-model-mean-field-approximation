#!/bin/bash


# Function to create a directory if it doesn't exist and add a .gitkeep file in order to allow git to tracking
create_directory() {
    if [ ! -d "$1" ]; then
        mkdir -p "$1"
        touch "$1/.gitkeep"
        echo "Created directory: $1"
    else
        echo "Directory already exists: $1"
    fi
}


# Define the list of directories~to be created. This must match the Makefile
directories=("obj" "obj/ini_file_parser")


# Create necessary directories
echo "Setting up project directories..."

# Loop through each element
for dir in "${directories[@]}"; do
  create_directory $dir
done

echo "Directories created successfully."
