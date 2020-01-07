#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TP63.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TP63 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TP73 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TP73-AS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPBG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPBGL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPCN1 
