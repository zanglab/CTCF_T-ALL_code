#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HIST1H3A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H3C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H3D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H3E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H3F 
