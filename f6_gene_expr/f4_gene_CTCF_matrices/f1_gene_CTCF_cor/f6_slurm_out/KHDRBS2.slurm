#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KHDRBS2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KHDRBS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KHDRBS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KHK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KHNYN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KHSRP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA0040 
