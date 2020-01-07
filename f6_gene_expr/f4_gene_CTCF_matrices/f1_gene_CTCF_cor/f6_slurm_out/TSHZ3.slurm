#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TSHZ3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TSHZ3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSKS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSKU 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSLP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSNARE1 
