#!/bin/bash
# Compile the E-QFT Mercury Standard Model paper

echo "Compiling E-QFT Mercury Standard Model paper..."

# Navigate to paper directory
cd ../paper

# Run pdflatex twice to resolve references
pdflatex E-QFT_Mercury_Standard_Model.tex
pdflatex E-QFT_Mercury_Standard_Model.tex

echo "Compilation complete. Output file: paper/E-QFT_Mercury_Standard_Model.pdf"