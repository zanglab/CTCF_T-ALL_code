#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o USP19.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i USP19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USP24 
