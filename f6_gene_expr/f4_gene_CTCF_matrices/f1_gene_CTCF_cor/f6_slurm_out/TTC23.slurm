#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TTC23.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TTC23 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTC23L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTC24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTC25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTC26 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTC27 
