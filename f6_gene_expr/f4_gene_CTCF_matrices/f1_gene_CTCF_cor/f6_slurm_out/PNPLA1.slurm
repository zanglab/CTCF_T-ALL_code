#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PNPLA1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PNPLA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PNPLA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PNPLA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PNPLA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PNPLA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PNPLA6 
