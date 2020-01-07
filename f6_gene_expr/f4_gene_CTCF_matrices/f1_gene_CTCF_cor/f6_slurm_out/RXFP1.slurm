#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RXFP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RXFP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RXFP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RXFP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RXFP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RXRA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RXRB 
