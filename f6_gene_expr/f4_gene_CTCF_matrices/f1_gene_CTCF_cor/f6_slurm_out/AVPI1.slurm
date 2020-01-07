#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AVPI1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AVPI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AVPR1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AVPR1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AVPR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AWAT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AWAT2 
