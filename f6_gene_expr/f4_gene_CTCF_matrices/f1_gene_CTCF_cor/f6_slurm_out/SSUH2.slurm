#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SSUH2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SSUH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SSX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SSX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SSX2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SSX2IP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SSX3 
