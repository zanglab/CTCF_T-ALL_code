#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HTATSF1P2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HTATSF1P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HTN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HTN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HTR1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HTR1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HTR1D 
