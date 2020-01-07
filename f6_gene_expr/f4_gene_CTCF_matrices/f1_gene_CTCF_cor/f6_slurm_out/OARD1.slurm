#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OARD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OARD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OAS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OAS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OAS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OASL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OAT 
