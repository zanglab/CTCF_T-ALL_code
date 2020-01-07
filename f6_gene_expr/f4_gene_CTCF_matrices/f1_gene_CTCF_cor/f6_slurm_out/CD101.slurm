#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CD101.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CD101 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD109 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD151 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD160 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD163 
