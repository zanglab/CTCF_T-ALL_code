#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o VCPKMT.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i VCPKMT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VCX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VCX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VCX3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VCX3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VCY1B 
