#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SEC14L6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SEC14L6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEC16A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEC16B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEC1P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEC22A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEC22B 
