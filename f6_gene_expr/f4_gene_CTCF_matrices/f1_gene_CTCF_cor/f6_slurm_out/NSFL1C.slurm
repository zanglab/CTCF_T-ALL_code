#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NSFL1C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NSFL1C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NSFP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NSG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NSL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NSMAF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NSMCE1 
