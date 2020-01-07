#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IMPA1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IMPA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IMPA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IMPACT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IMPAD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IMPDH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IMPDH2 
