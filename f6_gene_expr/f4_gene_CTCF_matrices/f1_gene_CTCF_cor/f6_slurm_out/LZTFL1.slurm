#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LZTFL1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LZTFL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LZTR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LZTS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LZTS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LZTS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i M1AP 
