#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LTBP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LTBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LTBP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LTBP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LTBP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LTBR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LTC4S 
