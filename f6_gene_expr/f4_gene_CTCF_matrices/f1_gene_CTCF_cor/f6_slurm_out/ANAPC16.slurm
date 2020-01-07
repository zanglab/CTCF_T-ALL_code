#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ANAPC16.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ANAPC16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANAPC1P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANAPC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANAPC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANAPC5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANAPC7 
