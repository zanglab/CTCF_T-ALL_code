#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RMND5A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RMND5A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RMND5B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNASE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNASE10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNASE11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RNASE13 
