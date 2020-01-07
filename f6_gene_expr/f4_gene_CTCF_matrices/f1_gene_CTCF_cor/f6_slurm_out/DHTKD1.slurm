#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DHTKD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DHTKD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DHX15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DHX16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DHX29 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DHX30 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DHX32 
