#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TERF2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TERF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TERF2IP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TERT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TES 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TESC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TESK1 
