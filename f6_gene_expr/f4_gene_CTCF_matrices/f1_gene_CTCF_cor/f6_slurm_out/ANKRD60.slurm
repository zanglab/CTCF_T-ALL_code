#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ANKRD60.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD60 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD61 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD62 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD63 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD65 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD66 
