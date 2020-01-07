#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IGF2BP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IGF2BP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGF2BP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGF2R 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGFALS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGFBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGFBP2 
