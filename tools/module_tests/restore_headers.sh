#!/bin/bash


# This script restores a file to its original state by replacing it with a previously created backup. 
# It is designed to work alongside scripts that temporarily modify files, ensuring that changes 
# can be safely reverted after they are no longer needed.
#
# The script takes a single argument: the name of the file to restore. It checks for the presence 
# of a backup file with the same name and a `.bak` extension in the current directory. If the backup 
# is found, the script replaces the modified file with the backup. If the backup is missing, it exits 
# with an error message.
#
# Usage:
# ./script_name.sh <file>
# Example: ./restore_headers.sh example.cpp


# Validate the input arguments
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <file>"
    exit 1
fi

# Assign the input argument to a variable
file="$1"
backup_file="${file}.bak"

# Ensure the backup file exists in the current directory
if [ ! -f "$backup_file" ]; then
    echo "Error: Backup file '${backup_file}' not found in the current directory."
    exit 1
fi

# Restore the original file
mv "$backup_file" "$file"

# Notify the user
echo "The file '${file}' has been restored to its original state."
