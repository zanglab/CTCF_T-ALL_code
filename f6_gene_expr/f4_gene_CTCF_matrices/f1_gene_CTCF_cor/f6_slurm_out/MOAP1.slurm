#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MOAP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MOAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MOB1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MOB1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MOB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MOB3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MOB3B 
