#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GSTA7P.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GSTA7P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GSTCD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GSTK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GSTM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GSTM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GSTM2P1 
