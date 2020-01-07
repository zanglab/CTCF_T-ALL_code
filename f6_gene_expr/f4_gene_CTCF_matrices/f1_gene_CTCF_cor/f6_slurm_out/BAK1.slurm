#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BAK1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BAK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BAMBI 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BANF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BANF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BANK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BANP 
