#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CXCL1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CXCL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CXCL10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CXCL11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CXCL12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CXCL13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CXCL14 
