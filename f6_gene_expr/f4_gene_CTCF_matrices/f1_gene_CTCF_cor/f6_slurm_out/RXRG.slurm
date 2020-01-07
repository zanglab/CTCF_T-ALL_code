#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RXRG.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RXRG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RYBP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RYK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RYR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RYR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RYR3 
