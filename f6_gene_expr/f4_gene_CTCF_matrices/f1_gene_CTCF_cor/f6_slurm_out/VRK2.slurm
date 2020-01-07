#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o VRK2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i VRK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VRK3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VRTN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VSIG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VSIG10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VSIG10L 
