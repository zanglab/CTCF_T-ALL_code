#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ANKRD50.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD50 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD52 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD53 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD54 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD55 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD6 
