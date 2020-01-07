#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RAB6C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB6C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB6C-AS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB7A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB7B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB8A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB8B 
