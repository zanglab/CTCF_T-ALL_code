#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ADGRL4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ADGRL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADGRV1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADH1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADH1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADH1C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADH4 
