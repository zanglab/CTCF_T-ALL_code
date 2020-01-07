#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ST8SIA2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ST8SIA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ST8SIA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ST8SIA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ST8SIA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ST8SIA6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STAB1 
