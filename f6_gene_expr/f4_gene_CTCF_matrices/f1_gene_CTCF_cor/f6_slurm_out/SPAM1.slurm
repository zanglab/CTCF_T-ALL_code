#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SPAM1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SPAM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPANXA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPANXB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPANXC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPANXD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPANXN1 
