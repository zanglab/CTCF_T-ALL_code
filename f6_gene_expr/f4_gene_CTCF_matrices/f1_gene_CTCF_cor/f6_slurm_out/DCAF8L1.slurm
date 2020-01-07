#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DCAF8L1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DCAF8L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCAF8L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCAKD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCANP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCBLD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCBLD2 
