#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LINC00993.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LINC00993 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LINC01085 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LINC01225 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LINC01347 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LINC01399 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LINC01529 
