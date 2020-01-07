#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CR1L.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CR1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRABP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRABP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRACR2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRACR2B 
