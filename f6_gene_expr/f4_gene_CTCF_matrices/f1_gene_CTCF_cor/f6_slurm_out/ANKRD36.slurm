#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ANKRD36.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD36 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD36B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD36BP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD36C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD37 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD39 
