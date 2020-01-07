#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o VBP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i VBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VCAM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VCAN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VCL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VCP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VCPIP1 
