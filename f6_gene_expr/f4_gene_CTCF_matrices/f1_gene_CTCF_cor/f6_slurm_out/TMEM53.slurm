#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMEM53.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM53 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM54 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM56 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM56-RWDD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM59 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM59L 
