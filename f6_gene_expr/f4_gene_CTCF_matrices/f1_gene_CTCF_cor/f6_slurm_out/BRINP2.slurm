#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BRINP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BRINP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BRINP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BRIP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BRIX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BRK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BRMS1 
