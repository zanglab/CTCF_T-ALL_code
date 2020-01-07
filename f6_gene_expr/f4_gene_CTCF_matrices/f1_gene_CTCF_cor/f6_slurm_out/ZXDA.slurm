#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZXDA.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZXDA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZXDB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZXDC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZYG11A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZYG11B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZYX 
