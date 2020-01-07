#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MOCS3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MOCS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MOG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MOGAT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MOGAT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MOGAT3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MOGS 
