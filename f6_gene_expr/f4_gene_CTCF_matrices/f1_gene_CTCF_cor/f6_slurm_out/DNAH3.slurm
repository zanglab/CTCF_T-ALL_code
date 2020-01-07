#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DNAH3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAH3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAH5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAH6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAH7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAH8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAH9 
