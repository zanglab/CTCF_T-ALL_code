#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NLRP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NLRP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NLRP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NLRP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NLRP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NLRP6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NLRP7 
