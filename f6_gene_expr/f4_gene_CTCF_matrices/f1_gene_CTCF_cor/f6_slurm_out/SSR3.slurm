#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SSR3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SSR3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SSR4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SSR4P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SSRP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SSSCA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SST 
