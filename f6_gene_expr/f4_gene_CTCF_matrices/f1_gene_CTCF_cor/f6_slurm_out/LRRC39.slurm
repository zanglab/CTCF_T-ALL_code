#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LRRC39.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC39 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC3C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC40 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRC41 
