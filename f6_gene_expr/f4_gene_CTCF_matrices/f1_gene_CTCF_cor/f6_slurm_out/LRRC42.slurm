#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LRRC42.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC42 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC43 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC45 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC46 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC47 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC49 
