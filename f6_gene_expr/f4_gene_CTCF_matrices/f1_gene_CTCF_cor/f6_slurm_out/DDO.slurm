#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DDO.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DDO 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDOST 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDRGK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDT 
