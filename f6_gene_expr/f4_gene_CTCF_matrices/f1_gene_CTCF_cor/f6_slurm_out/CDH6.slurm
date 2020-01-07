#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CDH6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CDH6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDH7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDH8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDH9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDHR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDHR2 
