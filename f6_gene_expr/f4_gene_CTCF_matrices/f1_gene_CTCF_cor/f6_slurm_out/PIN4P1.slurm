#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PIN4P1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PIN4P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PINK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PINLYP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PINX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIP4K2A 
