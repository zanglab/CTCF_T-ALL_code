#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RAB4B-EGLN2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB4B-EGLN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB5A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB5B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB5C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB6A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB6B 
