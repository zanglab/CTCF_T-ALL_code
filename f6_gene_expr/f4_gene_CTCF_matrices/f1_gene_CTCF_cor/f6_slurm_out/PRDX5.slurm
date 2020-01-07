#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PRDX5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PRDX5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRDX6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PREB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRELID1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRELID2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRELP 
