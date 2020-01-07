#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PRR23C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PRR23C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRR23D2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRR25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRR26 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRR27 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRR29 
