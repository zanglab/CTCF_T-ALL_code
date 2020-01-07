#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMC8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMC8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMCC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMCC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMCC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMCO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMCO2 
