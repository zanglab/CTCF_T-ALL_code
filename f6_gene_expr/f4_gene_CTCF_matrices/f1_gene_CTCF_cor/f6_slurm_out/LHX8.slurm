#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LHX8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LHX8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LHX9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LIAS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LIF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LIFR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LIG1 
