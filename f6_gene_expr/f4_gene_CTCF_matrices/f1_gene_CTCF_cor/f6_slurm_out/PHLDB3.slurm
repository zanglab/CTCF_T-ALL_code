#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PHLDB3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PHLDB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHLPP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHLPP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHOSPHO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHOSPHO2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHOX2A 
