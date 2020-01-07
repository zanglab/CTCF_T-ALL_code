#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNRF2P2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNRF2P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNRF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNRF4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZP3 
