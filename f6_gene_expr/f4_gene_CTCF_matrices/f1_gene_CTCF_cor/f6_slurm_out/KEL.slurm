#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KEL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KEL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KERA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KHDC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KHDC1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KHDC3L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KHDRBS1 
