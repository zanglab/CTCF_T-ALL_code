#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RCAN1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RCAN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RCAN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RCAN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RCBTB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RCBTB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RCC1 
