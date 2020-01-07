#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMEM265.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM265 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM30A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM30B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM31 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM33 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM37 
