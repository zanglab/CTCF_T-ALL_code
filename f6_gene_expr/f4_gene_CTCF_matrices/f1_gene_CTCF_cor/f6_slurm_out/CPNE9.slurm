#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CPNE9.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CPNE9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CPO 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CPOX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CPPED1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CPQ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CPS1 
