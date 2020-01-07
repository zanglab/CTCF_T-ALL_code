#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CSE1L.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CSE1L 
