#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CPLX4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CPLX4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CPM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CPN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CPN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CPNE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CPNE2 
