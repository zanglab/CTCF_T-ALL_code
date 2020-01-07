#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PDCD2L.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PDCD2L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDCD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDCD5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDCD6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDCD6IP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDCD6IPP2 
