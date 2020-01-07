#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LRRC34.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC34 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC36 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC37A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC37A11P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC37A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC37A3 
