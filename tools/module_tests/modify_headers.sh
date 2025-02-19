#!/bin/bash


# This script automates the process of modifying a specific text pattern in a file, 
# typically used to adjust include statements in C++ source or header files. 
# It takes three arguments: the file to modify, the text to search for, and the replacement text.
# 
# The script creates a backup of the original file to ensure changes can be reverted if needed. 
# It then applies the specified modifications using the `sed` command. 
# This tool is particularly useful in build workflows that require temporary changes 
# to file contents for testing or compatibility purposes.
# 
# Usage: 
# ./script_name.sh <file> <search> <replace>
# Example: ./modify_headers.sh example.cpp '#include "old.h"' '#include "new.h"'


# Validate the input arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <file> <search> <replace>"
    exit 1
fi

# Assign input arguments to variables
file="$1"
search="$2"
replace="$3"

# Backup the original file
backup_file="${file}.bak"
cp "$file" "$backup_file"

# Describe what the script is doing
echo "Executing a script to modify the include section in the file ${file}."
echo "Replacing '${search}' with '${replace}'."

# Modify the include statement
sed -i "s|$search|$replace|" "$file"

echo "Modification complete. Backup saved as ${backup_file}."
