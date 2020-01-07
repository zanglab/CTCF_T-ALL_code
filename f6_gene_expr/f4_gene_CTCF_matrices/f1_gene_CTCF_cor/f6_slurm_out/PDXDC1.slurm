#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PDXDC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PDXDC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDXDC2P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDXK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDXP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDYN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDZD11 
