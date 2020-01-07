#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CPVL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CPVL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CPXCR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CPXM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CPXM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CPZ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CR1 
