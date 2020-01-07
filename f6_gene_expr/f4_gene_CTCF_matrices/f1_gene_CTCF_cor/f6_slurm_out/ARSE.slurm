#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ARSE.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ARSE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARSF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARSG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARSH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARSI 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARSJ 
