#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AQP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AQP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AQP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AQP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AQP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AQP6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AQP7 
