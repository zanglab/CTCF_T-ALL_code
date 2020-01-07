#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LRRC2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC23 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC26 
