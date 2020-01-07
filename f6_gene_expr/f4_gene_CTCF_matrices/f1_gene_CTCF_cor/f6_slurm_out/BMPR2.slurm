#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BMPR2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BMPR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BMS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BMS1P17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BMS1P21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BMS1P4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BMX 
