#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RTKN2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RTKN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RTL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RTN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RTN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RTN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RTN4 
