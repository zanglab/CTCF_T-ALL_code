#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GDI2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GDI2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GDNF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GDPD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GDPD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GDPD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GDPD4 
