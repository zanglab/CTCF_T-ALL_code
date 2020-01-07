#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PRR18.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PRR18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRR19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRR20A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRR22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRR23A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRR23B 
