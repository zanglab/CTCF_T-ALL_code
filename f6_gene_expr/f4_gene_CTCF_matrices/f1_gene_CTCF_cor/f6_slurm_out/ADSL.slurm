#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ADSL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ADSL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADSS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADSSL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADTRP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AEBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AEBP2 
