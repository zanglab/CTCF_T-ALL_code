#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SRD5A3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SRD5A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SREBF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SREBF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SREK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SREK1IP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRF 
