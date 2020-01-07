#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PRKCG.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PRKCG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRKCH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRKCI 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRKCQ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRKCSH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRKCZ 
