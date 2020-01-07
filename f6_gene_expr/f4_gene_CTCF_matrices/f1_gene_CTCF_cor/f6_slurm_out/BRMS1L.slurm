#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BRMS1L.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BRMS1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BROX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BRPF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BRPF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BRS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BRSK1 
