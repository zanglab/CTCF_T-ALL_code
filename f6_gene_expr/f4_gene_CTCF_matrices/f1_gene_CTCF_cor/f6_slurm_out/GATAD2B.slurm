#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GATAD2B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GATAD2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GATB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GATC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GATM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GBA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GBA2 
