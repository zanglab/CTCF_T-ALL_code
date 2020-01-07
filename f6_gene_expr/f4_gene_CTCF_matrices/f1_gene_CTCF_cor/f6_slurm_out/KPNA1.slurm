#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KPNA1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KPNA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KPNA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KPNA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KPNA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KPNA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KPNA6 
