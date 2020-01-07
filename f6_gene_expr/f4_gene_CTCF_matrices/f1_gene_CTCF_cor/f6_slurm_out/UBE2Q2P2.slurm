#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o UBE2Q2P2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2Q2P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2QL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2R2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2S 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2T 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBE2U 
