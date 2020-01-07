#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LGALS7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LGALS7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LGALS7B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LGALS8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LGALS9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LGALS9B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LGALS9C 
