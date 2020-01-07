#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MARC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MARC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MARC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MARCH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MARCH10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MARCH11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MARCH2 
