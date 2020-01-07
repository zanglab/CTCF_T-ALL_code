#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AMELY.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AMELY 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AMER1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AMER2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AMER3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AMFR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AMH 
