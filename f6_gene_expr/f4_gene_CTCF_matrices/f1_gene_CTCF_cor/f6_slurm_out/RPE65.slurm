#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RPE65.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RPE65 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPEL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPGR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPGRIP1 
