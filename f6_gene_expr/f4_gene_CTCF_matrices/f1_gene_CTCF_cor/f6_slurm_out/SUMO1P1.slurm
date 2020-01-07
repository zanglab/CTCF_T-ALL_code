#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SUMO1P1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SUMO1P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SUMO1P3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SUMO2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SUMO3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SUMO4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SUN1 
