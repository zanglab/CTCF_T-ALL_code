#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SEPSECS.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SEPSECS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEPT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEPT10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEPT11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEPT12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEPT14 
