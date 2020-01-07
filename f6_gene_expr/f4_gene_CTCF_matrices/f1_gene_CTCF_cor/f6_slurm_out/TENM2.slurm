#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TENM2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TENM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TENM3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TENM4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TEP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TEPP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TERF1 
