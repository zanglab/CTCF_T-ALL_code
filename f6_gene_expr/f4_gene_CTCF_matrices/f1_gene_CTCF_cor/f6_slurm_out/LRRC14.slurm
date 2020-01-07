#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LRRC14.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC14B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC19 
