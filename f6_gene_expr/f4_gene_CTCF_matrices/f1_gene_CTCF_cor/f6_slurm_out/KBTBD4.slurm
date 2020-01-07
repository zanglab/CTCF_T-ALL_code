#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KBTBD4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KBTBD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KBTBD6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KBTBD7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KBTBD8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCMF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNA1 
