#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DPT.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DPT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPY19L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPY19L1P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPY19L1P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPY19L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPY19L2P1 
