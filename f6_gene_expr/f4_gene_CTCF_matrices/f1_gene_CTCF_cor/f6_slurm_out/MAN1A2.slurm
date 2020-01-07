#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MAN1A2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MAN1A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAN1B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAN1C1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAN2A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAN2A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAN2B1 
