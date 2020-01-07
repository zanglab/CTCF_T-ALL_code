#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ANKRD12.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD13A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD13B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD13C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD13D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD16 
