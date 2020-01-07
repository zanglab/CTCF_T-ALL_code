#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CDCA7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CDCA7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDCA7L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDCA8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDCP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDCP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDH1 
