#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o JAK2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i JAK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i JAK3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i JAKMIP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i JAKMIP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i JAKMIP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i JAM2 
