#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RETSAT.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RETSAT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i REV1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i REV3L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i REXO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i REXO2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i REXO4 
