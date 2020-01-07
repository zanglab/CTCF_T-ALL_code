#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MAGEA2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGEA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGEA2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGEA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGEA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGEA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGEA6 
