#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RERGL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RERGL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RESP18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i REST 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RET 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RETN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RETNLB 
