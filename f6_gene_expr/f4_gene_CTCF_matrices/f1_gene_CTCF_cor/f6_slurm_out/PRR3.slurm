#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PRR3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PRR3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRR30 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRR32 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRR34 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRR35 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRR36 
