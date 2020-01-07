#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o VPS37A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i VPS37A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VPS37B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VPS37C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VPS37D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VPS39 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VPS41 
