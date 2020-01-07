#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMEM220.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM220 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM221 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM222 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM223 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM225 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM229A 
