#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o JUNB.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i JUNB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i JUND 
time python ../7_cor_CTCF_binding_GeneExpr.py -i JUP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KAAG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KALRN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KANK1 
