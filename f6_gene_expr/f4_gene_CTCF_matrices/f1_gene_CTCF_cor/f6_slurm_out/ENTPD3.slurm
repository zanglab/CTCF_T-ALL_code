#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ENTPD3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ENTPD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ENTPD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ENTPD5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ENTPD6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ENTPD7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ENTPD8 
