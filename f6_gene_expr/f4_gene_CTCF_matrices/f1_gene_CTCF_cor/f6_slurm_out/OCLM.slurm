#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OCLM.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OCLM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OCLN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OCM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OCM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OCRL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OCSTAMP 
