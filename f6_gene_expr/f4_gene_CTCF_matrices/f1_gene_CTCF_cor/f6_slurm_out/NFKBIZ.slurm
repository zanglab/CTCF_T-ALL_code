#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NFKBIZ.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NFKBIZ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NFRKB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NFS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NFU1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NFX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NFXL1 
