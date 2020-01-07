#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RAB3A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB3C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB3D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB3GAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB3GAP2 
