#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MYBPC2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MYBPC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYBPC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYBPH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYBPHL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYCBP 
