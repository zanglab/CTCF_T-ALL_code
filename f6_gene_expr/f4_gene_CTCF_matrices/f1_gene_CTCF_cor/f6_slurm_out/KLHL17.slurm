#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KLHL17.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLHL22 
