#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PRKACG.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PRKACG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRKAG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRKAG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRKAG3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRKAR1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRKAR1B 
