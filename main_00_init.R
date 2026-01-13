### ======================================================================== ###
rm(list=ls())
set.seed(12345)
### ======================================================================== ###


### ======================================================================== ###
### check and install required packages ###
### ======================================================================== ###
required_packages <- c("readr", "xts", "Rsolnp", "matrixStats",  
                       "expm", "matrixcalc")
new_packages <- required_packages[!required_packages %in% installed.packages()]
new_packages
if(length(new_packages)) install.packages(new_packages)
### ======================================================================== ###


### ======================================================================== ###
### load required packages ###
### ======================================================================== ###
library("readr")
library("xts")
library("Rsolnp")
library("matrixStats")
library("expm")
library("matrixcalc")
### ======================================================================== ###