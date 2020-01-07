#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FKBP1AP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FKBP1AP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FKBP1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FKBP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FKBP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FKBP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FKBP5 
