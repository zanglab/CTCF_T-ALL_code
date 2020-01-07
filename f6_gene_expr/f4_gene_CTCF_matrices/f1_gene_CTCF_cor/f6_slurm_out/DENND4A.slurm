#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DENND4A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DENND4A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DENND4B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DENND4C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DENND5A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DENND5B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DENND6A 
