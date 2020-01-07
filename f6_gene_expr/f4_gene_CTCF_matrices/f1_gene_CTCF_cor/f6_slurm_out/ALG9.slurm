#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ALG9.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ALG9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALKBH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALKBH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALKBH3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALKBH4 
