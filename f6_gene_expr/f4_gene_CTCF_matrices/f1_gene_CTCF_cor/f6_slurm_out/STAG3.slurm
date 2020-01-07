#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o STAG3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i STAG3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STAG3L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STAG3L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STAG3L3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STAG3L4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STAG3L5P 
