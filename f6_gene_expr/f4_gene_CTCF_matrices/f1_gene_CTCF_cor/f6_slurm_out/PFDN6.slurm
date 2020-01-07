#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PFDN6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PFDN6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PFKFB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PFKFB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PFKFB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PFKFB4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PFKL 
