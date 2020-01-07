#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GLRB.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GLRB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLRX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLRX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLRX3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLRX5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLS 
