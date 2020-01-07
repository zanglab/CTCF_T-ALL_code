#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o STYK1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i STYK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STYX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STYXL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SUB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SUCLA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SUCLG1 
