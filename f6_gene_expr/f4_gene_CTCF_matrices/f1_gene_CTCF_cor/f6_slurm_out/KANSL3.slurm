#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KANSL3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KANSL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KANTR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KARS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KAT2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KAT2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KAT5 
