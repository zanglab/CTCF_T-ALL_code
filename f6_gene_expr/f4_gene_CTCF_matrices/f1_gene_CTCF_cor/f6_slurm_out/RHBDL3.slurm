#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RHBDL3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RHBDL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RHBG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RHCE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RHCG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RHD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RHEB 
