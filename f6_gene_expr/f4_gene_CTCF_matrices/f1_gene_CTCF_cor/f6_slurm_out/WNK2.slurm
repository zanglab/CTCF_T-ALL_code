#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o WNK2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i WNK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WNK3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WNK4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WNT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WNT10A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WNT10B 
