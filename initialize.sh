#!/bin/bash

echo "Conda environment installation directory: $1";
echo "Conda source path: $2"

if [ ! -d "$1/envs" ]; then
    mkdir $1/envs
fi
for file in "$(pwd)/environment_configs"/*; do
    # Check if the current item is a file
    if [ -f "$file" ]; then
        # Create a conda environment for every file in the environment_configs directory if it does not already exist
        env_name="${file##*/}"
        env_name="${env_name%.yml}"
        # Check if the environment already exists
        if [ ! -d "$1/envs/$env_name" ]; then
            echo "Creating environment using : $file"
            CONDA_PKGS_DIRS=$(mktemp -d) mamba env create --prefix="$1/envs/$env_name" --file="$file"
            if [[ $file == *'boostdiff'* ]]; then
                # special installation steps needed for boostdiff
                current_directory=$(pwd)
                source $2/activate "$1/envs/$env_name"
                cd "$(pwd)/boostdiff_inference"
                pip install .
                cd $current_directory
                source $2/deactivate
            fi
            if [[ $file == *'shiny'* ]]; then
                # special installation steps needed for boostdiff
                current_directory=$(pwd)
                source $2/activate "$1/envs/$env_name"
                Rscript install_shiny_packages.R
                source $2/deactivate
            fi
            # Add your own special installation steps for an environment if needed:
        fi
    fi
done
