#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LILRB5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LILRB5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LILRP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LIM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LIMA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LIMCH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LIMD1 
