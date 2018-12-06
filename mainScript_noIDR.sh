#!/bin/bash 
# 
#$ -cwd 
#$ -V 
#$ -j y 
#$ -S /bin/bash 
# 
 
R --vanilla <./code/mainScript_noIDR.R> mainScript_noIDR.Rout
