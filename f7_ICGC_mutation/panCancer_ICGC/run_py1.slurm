#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p largemem 
#SBATCH -A zanglab
#SBATCH -o out_py1.log

#Run program
time python py1_mutation_around_CTCFs.py

