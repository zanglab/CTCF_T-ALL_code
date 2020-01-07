#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FRMPD2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FRMPD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FRMPD2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FRMPD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FRMPD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FRRS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FRRS1L 
