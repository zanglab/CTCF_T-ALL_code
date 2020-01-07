#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SERP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SERP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERPINA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERPINA10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERPINA11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERPINA12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERPINA13P 
