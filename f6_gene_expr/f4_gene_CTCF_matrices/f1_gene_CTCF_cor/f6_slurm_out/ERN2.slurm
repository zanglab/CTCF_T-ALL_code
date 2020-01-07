#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ERN2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ERN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERO1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERO1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERP27 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERP29 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERP44 
