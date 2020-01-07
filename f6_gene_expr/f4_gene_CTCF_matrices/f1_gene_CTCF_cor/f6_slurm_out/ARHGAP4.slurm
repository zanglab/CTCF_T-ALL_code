#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ARHGAP4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGAP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGAP40 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGAP42 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGAP44 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGAP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGAP6 
