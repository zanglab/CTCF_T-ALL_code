#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMEM38A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM38A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM38B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM39A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM39B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM40 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM41A 
