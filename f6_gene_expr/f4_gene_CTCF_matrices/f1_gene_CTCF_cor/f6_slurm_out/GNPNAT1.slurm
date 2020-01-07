#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GNPNAT1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GNPNAT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNPTAB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNPTG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNRH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNRH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNRHR 
