#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LRRC8C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC8C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC8D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC8E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRCC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRD1 
