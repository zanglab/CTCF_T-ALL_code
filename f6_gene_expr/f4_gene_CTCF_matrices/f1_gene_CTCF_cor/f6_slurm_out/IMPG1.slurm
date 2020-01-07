#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IMPG1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IMPG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IMPG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INAFM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INAFM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INCA1 
