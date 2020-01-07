#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SS18L2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SS18L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SSB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SSBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SSBP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SSBP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SSBP4 
