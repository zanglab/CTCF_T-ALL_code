#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR10C1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR10C1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR10G2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR10G3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR10G4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR10G7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR10G8 
