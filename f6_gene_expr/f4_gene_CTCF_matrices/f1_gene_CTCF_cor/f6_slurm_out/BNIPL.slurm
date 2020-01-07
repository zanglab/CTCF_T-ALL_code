#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BNIPL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BNIPL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BOC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BOD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BOD1L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BOD1L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BOK 
