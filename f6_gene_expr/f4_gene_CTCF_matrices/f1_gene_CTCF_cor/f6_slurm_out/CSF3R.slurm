#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CSF3R.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CSF3R 
