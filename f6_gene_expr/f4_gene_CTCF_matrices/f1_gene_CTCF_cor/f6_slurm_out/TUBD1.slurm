#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TUBD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TUBD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TUBE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TUBG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TUBG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TUBGCP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TUBGCP3 
