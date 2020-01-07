#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LHX1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LHX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LHX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LHX3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LHX4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LHX5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LHX6 
