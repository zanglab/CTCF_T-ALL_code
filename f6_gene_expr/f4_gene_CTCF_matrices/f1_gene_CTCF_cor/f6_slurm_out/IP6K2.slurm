#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IP6K2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IP6K2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IP6K3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IPCEF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IPMK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IPO11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IPO13 
