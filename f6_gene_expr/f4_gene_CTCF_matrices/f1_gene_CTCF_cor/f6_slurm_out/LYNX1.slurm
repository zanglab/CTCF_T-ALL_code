#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LYNX1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LYNX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LYPD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LYPD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LYPD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LYPD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LYPD5 
