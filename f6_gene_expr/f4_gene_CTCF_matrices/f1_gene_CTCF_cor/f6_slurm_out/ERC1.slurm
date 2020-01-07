#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ERC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ERC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERCC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERCC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERCC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERCC4 
