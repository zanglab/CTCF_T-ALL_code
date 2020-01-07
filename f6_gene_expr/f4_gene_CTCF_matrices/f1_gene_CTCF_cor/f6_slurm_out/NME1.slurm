#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NME1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NME1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NME1-NME2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NME2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NME3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NME4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NME5 
