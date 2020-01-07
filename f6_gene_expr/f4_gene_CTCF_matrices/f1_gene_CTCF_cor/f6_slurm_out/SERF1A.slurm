#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SERF1A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SERF1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERF1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERGEF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERHL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERHL2 
