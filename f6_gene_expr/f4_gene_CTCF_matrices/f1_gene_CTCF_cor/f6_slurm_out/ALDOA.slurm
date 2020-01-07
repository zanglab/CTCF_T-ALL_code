#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ALDOA.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ALDOA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALDOB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALDOC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALG10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALG10B 
