#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LRRC27.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC27 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC28 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC29 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC31 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC32 
