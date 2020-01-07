#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BCAS3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BCAS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BCAS4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BCAT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BCAT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BCCIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BCDIN3D 
