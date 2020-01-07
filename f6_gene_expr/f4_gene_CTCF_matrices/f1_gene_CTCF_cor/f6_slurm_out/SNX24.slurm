#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SNX24.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SNX24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNX25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNX27 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNX29 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNX29P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNX29P2 
