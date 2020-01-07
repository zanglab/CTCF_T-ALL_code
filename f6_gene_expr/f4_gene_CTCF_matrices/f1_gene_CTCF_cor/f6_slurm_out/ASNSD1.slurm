#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ASNSD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ASNSD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASPA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASPDH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASPG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASPH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASPHD1 
