#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o APTX.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i APTX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AQP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AQP10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AQP11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AQP12A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AQP12B 
