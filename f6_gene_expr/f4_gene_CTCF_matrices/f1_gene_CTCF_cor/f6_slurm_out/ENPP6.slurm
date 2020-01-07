#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ENPP6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ENPP6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ENPP7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ENSA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ENTHD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ENTPD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ENTPD2 
