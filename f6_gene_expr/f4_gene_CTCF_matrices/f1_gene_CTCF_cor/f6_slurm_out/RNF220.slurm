#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RNF220.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF220 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF222 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF223 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF224 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF225 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNF24 
