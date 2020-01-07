#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CTSG.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CTSG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTSH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTSK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTSL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTSL3P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTSLP2 
