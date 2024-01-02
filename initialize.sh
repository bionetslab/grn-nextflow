#!/bin/bash

echo "Conda environment directory: $1";
echo "Conda source path: $2"

for file in "$1/environment_configs"/*; do
    # Check if the current item is a file (not a directory)
    if [ -f "$file" ]; then
        env_name="${file##*/}"
        env_name="${env_name%.yml}"
        if [ ! -d "$1/envs/$env_name" ]; then
            echo "Creating environment using : $file"
            CONDA_PKGS_DIRS=$(mktemp -d) mamba env create --prefix="$1/envs/$env_name" --file="$file"
            if [[ $file == *'boostdiff'* ]]; then
                current_directory=$(pwd)
                source $2/activate "$1/envs/$env_name"
                cd "$(pwd)/boostdiff_inference"
                pip install .
                cd $current_directory
                source $2/deactivate
            fi
        fi
        # mamba env create --prefix="$1/envs/$env_name" --file="$file"
        
        # Add your processing logic here
        # For example, you can read the contents of the file, perform some operations, etc.
        
    fi
done