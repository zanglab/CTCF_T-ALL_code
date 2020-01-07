#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RARA.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RARA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RARB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RARG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RARRES1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RARRES2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RARRES3 
