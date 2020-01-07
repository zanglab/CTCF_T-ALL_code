#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o A4GALT.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i A4GALT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i A4GNT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AAAS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AACS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AACSP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AADAC 
