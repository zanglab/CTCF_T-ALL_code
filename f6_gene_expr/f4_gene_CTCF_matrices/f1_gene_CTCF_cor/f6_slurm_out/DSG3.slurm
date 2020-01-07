#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DSG3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DSG3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DSG4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DSN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DSP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DSPP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DST 
