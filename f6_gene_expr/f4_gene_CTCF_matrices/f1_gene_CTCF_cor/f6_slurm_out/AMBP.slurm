#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AMBP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AMBP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AMBRA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AMD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AMDHD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AMDHD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AMELX 
