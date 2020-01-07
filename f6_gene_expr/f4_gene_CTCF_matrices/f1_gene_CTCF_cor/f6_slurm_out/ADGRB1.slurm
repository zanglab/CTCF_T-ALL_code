#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ADGRB1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ADGRB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADGRB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADGRB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADGRD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADGRE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADGRE2 
