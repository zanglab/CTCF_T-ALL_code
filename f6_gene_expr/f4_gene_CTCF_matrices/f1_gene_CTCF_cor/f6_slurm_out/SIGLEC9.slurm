#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SIGLEC9.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SIGLEC9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIGLECL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIGMAR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIK3 
