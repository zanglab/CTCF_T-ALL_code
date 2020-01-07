#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o APBB1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i APBB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APBB1IP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APBB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APBB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APC2 
