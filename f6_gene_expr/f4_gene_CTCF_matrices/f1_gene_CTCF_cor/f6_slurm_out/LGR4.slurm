#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LGR4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LGR4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LGR5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LGR6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LGSN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LHB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LHCGR 
