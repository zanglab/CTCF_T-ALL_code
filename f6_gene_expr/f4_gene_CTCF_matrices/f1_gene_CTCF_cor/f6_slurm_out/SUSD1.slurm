#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SUSD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SUSD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SUSD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SUSD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SUSD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SUSD5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SUSD6 
