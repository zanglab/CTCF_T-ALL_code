#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HVCN1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HVCN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HYAL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HYAL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HYAL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HYAL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HYALP1 
