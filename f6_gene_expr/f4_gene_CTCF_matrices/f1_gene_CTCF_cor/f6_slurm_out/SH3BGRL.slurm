#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SH3BGRL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SH3BGRL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SH3BGRL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SH3BGRL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SH3BP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SH3BP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SH3BP4 
