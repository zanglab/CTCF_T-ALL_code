#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ACSL4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ACSL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACSL5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACSL6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACSM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACSM2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACSM2B 
