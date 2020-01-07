#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TNFRSF11A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TNFRSF11A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNFRSF11B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNFRSF12A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNFRSF13B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNFRSF13C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNFRSF14 
