#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DLX1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DLX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DLX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DLX3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DLX4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DLX5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DLX6 
