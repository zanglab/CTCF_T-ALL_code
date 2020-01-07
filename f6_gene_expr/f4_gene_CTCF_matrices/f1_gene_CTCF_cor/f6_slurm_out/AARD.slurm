#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AARD.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AARD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AARS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AARS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AARSD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AASDH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AASDHPPT 
