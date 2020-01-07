#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMEM159.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM159 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM160 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM161A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM161B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM163 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM164 
