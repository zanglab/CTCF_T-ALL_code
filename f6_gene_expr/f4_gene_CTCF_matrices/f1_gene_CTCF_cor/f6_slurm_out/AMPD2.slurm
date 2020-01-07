#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AMPD2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AMPD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AMPD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AMPH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AMT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AMTN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AMY1A 
