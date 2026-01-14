#!/bin/bash

# Optional: Exit if any command fails
#set -e

# Define container and script path (modify as needed)
CONTAINER="/home/dillin01/cmr/cmrsim_jupyter.sif"
SCRIPT="/home/dillin01/cmr/simulation_batch.py"


PARENT_DIR="/home/dillin01/cmr/json"

for dir in "$PARENT_DIR"/*/; do
	dir_name=$(basename "$dir")
	echo "Found directory: $dir_name"

	#for rpm in 50 200 500 -50 -200 -500; do
	for rpm in -50; do
		echo "Running $rpm rpm case."

		singularity exec --nv "$CONTAINER" python "$SCRIPT" --seq "$dir_name" --dir "$PARENT_DIR" --rpm "$rpm"
	done
done
