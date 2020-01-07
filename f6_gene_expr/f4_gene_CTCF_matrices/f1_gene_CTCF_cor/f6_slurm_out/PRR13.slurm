#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PRR13.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PRR13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRR14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRR14L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRR15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRR15L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRR16 
