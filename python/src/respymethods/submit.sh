#!/bin/bash
#SBATCH --job-name=test-sim-respymethods
#SBATCH --output=%j.out
#SBATCH --time=48:00:00
#SBATCH --partition=zen5
#SBATCH --cpus-per-task=93
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elio.balestrieri@gmail.com

export OPENBLAS_NUM_THREADS=1

ml uv
uv run sim-sandbox.py


