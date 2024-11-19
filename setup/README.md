# Setup

## Setting Up Project Directories

The directory structure required for the project must exist before running the make command, as it relies on specific folders (e.g., obj and bin). Failure to create these directories beforehand will result in a build failure.

To automate this step, use the setup/create_directories.sh script, which ensures the necessary directories are created and adds a .gitkeep file to each one. The .gitkeep file allows Git to track these otherwise-empty directories, maintaining the static structure of the code for easier version control.
Usage

Run the script from the project root directory:
```
bash setup/create_directories.sh
```
This script:
- Creates the required directories if they do not exist.
- Adds a .gitkeep file to ensure the folders are tracked by Git.