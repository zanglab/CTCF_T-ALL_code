#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NPHP3-ACAD11.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NPHP3-ACAD11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPHP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPHS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPHS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPIPA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPIPA2 
