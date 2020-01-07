#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMEM176A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM176A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM176B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM177 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM178A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM178B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM179 
