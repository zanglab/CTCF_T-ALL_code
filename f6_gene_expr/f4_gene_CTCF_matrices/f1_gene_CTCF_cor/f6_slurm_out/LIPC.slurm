#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LIPC.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LIPC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LIPE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LIPF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LIPG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LIPH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LIPI 
