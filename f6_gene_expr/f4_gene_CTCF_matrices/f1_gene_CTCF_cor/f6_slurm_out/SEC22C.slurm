#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SEC22C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SEC22C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEC23A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEC23B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEC23IP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEC24A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEC24B 
