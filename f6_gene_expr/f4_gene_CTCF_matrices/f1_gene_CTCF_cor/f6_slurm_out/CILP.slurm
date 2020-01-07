#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CILP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CILP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CILP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CINP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CIPC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CIR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CIRBP 
