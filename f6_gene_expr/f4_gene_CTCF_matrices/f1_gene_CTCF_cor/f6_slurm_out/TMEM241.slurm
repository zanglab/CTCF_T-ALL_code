#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMEM241.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM241 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM242 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM243 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM244 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM245 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM246 
