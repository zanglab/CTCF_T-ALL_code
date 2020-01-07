#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DNAJC9.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJC9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNALI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNASE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNASE1L1 
