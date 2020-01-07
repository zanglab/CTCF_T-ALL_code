#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SDK2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SDK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SDR16C5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SDR16C6P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SDR39U1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SDR42E1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SDR9C7 
