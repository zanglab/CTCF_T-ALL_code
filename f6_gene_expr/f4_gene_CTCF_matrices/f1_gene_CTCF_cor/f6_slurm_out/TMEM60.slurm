#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMEM60.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM60 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM61 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM62 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM63A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM63B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM63C 
