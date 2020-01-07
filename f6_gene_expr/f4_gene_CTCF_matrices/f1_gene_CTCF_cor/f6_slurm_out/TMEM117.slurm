#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMEM117.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM117 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM119 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM120A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM120B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM121 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM123 
