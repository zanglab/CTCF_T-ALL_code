#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ANKRD33.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD33 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD33B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD34A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD34B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD34C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD35 
