#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BASP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BASP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BASP1P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BATF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BATF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BATF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BAX 
