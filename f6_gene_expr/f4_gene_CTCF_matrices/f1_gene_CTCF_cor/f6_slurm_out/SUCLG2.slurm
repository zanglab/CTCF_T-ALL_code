#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SUCLG2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SUCLG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SUCNR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SUCO 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SUDS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SUFU 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SUGCT 
