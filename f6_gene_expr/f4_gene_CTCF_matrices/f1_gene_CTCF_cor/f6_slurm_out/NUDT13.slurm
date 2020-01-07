#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NUDT13.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NUDT13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUDT14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUDT15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUDT16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUDT16L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUDT16P1 
