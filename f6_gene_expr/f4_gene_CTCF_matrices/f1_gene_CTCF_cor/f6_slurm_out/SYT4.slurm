#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SYT4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SYT4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYT5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYT6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYT7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYT8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYT9 
