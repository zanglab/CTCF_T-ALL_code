#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CTAG1A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CTAG1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTAG1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTAG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTAGE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTAGE10P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTAGE11P 
