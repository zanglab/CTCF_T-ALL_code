#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMC2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMC5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMC6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMC7 
