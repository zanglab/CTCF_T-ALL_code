#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PHOX2B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PHOX2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHPT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHRF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHTF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHTF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHYH 
