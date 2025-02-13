#!/bin/bash

# Exit script on error
set -e

# Change directory to script directory
HOSTILEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $HOSTILEDIR

# Define Conda environment name
ENV_NAME="hostile"

# Check if Conda is installed
if ! command -v conda &> /dev/null; then
    echo "Error: Conda is not installed. Please install Miniconda or Anaconda."
    exit 1
fi

# Create Conda environment if it doesn't exist
if ! conda info --envs | awk '{print $1}' | grep -Fxq "$ENV_NAME"; then
    echo "Creating Conda environment: $ENV_NAME"
    conda create -y -n "$ENV_NAME" python=3.10 hostile
else
    echo "Conda environment $ENV_NAME already exists."
fi

# Activate the environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

# Verify installation
if ! command -v hostile &> /dev/null; then
    echo "Error: Hostile is not installed correctly."
    exit 1
fi

# Set cache  directory for index to be saved
export HOSTILE_CACHE_DIR=${HOSTILEDIR}

# Fetch Hostile index
echo "Fetching Hostile index..."
hostile index fetch

conda deactivate

echo "✅ Hostile index fetch completed successfully!"
