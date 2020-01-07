#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RPL35A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL35A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL36 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL36A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL36A-HNRNPH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL36AL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL37 
