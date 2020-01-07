#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ANKRD40.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD40 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD42 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD44 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD45 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD46 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD49 
