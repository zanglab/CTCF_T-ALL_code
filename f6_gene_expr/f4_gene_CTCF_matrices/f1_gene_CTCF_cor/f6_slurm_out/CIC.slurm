#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CIC.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CIC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CIDEA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CIDEB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CIDEC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CIDECP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CIITA 
