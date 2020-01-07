#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ADAT3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAT3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADCK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADCK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADCK5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADCY1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADCY10 
