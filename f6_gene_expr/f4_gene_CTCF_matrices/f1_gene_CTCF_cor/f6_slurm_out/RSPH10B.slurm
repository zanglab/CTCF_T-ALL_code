#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RSPH10B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RSPH10B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RSPH10B2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RSPH14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RSPH3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RSPH4A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RSPH6A 
