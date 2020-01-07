#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GPN1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GPN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPNMB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR101 
