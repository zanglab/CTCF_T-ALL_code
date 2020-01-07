#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TOB1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TOB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TOB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TOB2P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TOE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TOLLIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TOM1 
