#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FKRP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FKRP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FKTN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FLAD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FLCN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FLG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FLG2 
