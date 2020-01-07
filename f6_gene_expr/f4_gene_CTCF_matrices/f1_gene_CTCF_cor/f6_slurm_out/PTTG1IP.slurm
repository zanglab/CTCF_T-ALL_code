#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PTTG1IP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PTTG1IP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTTG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTTG3P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTX3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTX4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PUF60 
