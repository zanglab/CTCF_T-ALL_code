#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RAD51.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RAD51 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAD51AP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAD51AP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAD51B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAD51C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAD51D 
