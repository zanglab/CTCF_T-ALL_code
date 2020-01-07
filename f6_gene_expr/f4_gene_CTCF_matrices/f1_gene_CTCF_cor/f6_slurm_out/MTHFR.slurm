#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MTHFR.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MTHFR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTHFS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTHFSD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTIF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTIF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTM1 
