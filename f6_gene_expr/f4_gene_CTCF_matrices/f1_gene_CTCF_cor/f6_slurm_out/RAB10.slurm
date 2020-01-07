#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RAB10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB11A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB11B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB11FIP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB11FIP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB11FIP3 
