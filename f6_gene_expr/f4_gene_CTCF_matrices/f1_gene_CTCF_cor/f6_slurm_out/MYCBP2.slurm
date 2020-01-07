#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MYCBP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MYCBP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYCBPAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYCL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYCN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYCT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYD88 
