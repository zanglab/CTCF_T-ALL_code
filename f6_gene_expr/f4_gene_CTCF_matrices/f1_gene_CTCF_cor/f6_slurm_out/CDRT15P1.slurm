#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CDRT15P1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CDRT15P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDRT4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDSN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDT1 
