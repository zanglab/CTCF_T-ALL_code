#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HTR3B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HTR3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HTR3C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HTR3D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HTR3E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HTR4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HTR5A 
