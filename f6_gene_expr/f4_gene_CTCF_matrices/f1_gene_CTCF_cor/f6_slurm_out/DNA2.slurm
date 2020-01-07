#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DNA2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DNA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAAF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAAF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAAF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAAF5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAH1 
