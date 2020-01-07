#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NECAP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NECAP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEDD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEDD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEDD4L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEDD8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEDD8-MDP1 
