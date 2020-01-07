#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CPEB3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CPEB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CPEB4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CPED1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CPLX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CPLX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CPLX3 
