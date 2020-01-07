#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RAB17.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB20 
