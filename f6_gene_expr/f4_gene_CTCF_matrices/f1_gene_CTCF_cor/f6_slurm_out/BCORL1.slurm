#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BCORL1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BCORL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BCORP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BCR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BCRP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BCS1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BDH1 
