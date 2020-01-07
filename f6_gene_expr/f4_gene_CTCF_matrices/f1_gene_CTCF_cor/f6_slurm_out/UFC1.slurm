#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o UFC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i UFC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UFL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UFM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UFSP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UFSP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UGCG 
