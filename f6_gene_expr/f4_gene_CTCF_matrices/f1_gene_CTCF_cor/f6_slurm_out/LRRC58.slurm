#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LRRC58.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC58 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC59 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC61 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC63 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC66 
