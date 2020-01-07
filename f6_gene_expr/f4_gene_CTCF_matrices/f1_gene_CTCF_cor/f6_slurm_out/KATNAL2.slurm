#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KATNAL2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KATNAL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KATNB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KATNBL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KATNBL1P6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KAZALD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KAZN 
