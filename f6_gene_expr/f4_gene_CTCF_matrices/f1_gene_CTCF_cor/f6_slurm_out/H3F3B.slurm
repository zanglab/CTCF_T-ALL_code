#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o H3F3B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i H3F3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i H3F3C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i H6PD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HAAO 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HABP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HABP4 
