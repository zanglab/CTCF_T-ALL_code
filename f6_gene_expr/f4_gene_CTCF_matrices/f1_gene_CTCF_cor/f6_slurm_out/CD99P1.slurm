#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CD99P1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CD99P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDADC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDAN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDC123 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDC14A 
