#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o VPREB1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i VPREB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VPREB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VPS11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VPS13A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VPS13B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VPS13C 
