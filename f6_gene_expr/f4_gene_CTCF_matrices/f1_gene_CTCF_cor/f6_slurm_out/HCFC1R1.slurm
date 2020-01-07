#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HCFC1R1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HCFC1R1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HCFC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HCG4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HCG4B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HCK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HCLS1 
