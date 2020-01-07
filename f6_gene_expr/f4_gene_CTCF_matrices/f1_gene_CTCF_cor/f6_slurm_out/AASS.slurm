#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AASS.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AASS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AATF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AATK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABAT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABCA10 
