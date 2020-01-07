#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HSPE1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HSPE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSPE1-MOB4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSPG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSPH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HTATIP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HTATSF1 
