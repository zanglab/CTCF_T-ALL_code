#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DPP4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DPP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPP6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPP7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPP8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPP9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPPA2 
