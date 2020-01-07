#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SKI.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SKI 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SKIDA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SKIL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SKIV2L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SKOR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SKOR2 
