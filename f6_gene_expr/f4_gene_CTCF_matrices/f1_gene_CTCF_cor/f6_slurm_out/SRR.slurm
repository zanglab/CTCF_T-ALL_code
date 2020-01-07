#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SRR.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SRR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRRD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRRM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRRM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRRM3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRRM4 
