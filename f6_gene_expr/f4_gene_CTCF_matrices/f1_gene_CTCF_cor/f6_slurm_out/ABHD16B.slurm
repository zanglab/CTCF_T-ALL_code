#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ABHD16B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ABHD16B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABHD17A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABHD17B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABHD17C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABHD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABHD3 
