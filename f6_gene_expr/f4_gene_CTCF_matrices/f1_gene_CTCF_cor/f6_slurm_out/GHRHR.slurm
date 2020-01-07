#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GHRHR.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GHRHR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GHRL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GHSR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GID4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GID8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GIF 
