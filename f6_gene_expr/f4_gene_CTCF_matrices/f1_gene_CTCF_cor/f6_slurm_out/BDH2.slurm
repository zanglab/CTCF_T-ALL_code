#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BDH2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BDH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BDKRB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BDKRB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BDNF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BDP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BEAN1 
