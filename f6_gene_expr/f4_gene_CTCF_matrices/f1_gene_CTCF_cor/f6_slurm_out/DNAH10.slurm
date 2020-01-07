#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DNAH10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAH10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAH11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAH12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAH14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAH17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAH2 
