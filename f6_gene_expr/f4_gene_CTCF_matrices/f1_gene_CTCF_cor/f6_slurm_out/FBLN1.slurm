#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FBLN1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FBLN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBLN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBLN5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBLN7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBN2 
