#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KIAA1614.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA1614 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA1671 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA1755 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA1841 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA1958 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA2012 
