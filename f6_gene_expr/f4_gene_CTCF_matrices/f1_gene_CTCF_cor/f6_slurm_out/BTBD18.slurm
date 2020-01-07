#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BTBD18.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BTBD18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BTBD19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BTBD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BTBD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BTBD6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BTBD7 
