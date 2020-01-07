#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NCAPD3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NCAPD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCAPG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCAPG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCAPH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCAPH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCBP1 
