#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMEM179B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM179B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM181 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM182 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM183A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM183B 
