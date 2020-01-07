#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TGDS.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TGDS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TGFA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TGFB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TGFB1I1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TGFB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TGFB3 
