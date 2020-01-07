#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DPY19L2P2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DPY19L2P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPY19L2P3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPY19L2P4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPY19L3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPY19L4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPY30 
