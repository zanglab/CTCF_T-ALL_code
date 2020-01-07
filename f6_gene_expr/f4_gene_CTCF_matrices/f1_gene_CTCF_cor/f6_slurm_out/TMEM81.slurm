#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMEM81.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM81 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM82 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM86A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM86B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM87A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM87B 
