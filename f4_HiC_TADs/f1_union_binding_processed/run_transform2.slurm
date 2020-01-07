#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 48:00:00
#SBATCH -p gpu 
#SBATCH -A zanglab
#SBATCH -o transform2.out

time python hicpro_matrix_to_jurcer_matrix_format.py -d PrEC
