#!/bin/bash
# Run Mercury E-QFT Simulation and Generate Figures
# This script runs the Mercury orbit simulation and generates figures for the paper

echo "========== E-QFT Mercury Orbit Simulation =========="
echo "Starting simulation process..."

# Make sure scripts are executable
chmod +x ../src/mercury_eqft_simulation.py
chmod +x ../src/generate_figures.py

# Navigate to source directory
cd ../src

# Run the main simulation
echo "Running main Mercury orbit simulation..."
python3 mercury_eqft_simulation.py

# Generate figures for the paper
echo "Generating figures for the paper..."
python3 generate_figures.py

# Copy the latest generated figures to the figures directory
echo "Copying figures to paper_assets directory..."
cp mercury_orbit_eqft.png ../figures/
cp figures_*/*.png ../figures/
cp mercury_precession_table.tex ../paper/paper_assets/

echo "Done! All simulation assets have been generated."
echo "Results can be found in the 'figures' directory."