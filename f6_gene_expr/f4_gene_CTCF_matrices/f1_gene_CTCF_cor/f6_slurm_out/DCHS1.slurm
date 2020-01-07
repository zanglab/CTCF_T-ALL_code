#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DCHS1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DCHS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCHS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCLK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCLK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCLK3 
