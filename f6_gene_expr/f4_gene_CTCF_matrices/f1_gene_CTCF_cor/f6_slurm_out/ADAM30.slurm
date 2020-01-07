#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ADAM30.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAM30 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAM32 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAM33 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAM3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAM5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAM6 
