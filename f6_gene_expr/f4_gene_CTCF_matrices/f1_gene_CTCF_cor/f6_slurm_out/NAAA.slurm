#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NAAA.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NAAA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAALAD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAALADL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAALADL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAB2 
