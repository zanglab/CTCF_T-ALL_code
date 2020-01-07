#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ADAMTS18.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAMTS18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAMTS19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAMTS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAMTS20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAMTS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAMTS4 
