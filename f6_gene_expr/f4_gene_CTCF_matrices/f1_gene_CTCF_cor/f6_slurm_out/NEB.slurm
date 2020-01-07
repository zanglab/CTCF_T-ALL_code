#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NEB.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NEB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEBL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NECAB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NECAB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NECAB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NECAP1 
