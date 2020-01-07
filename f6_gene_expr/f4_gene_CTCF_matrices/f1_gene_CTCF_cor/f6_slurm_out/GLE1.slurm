#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GLE1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GLE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLI2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLI3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLI4 
