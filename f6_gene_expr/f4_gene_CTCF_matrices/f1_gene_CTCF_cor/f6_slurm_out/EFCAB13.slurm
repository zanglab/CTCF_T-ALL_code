#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EFCAB13.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EFCAB13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EFCAB14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EFCAB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EFCAB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EFCAB5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EFCAB6 
