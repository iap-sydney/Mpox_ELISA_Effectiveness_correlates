# -------------------------------------------------------------------------
# Load all packages necessary to run the analaysis


library(rstan)
library(bayesplot)
library(tidyverse)
library(loo)
library(ggpubr)
library(glue)
library(latex2exp)
library(nleqslv)


# Global Parameters -------------------------------------------------------
source('R_functions/decay_titer.R')

#Specify Stan sampling parameters
warmup<-2000
sampling<-7000


# Define Colors used in figures -------------------------------------------
colD<-c("dodgerblue2","mediumpurple3","darkorange2")
colD2<-c("dodgerblue2","mediumpurple3","forestgreen","darkorange2")
colV<-c("darkorange2","dodgerblue2","mediumpurple3")
col_scale<-c('brown3',
             'blue3',
             'aquamarine3',
             'chartreuse3',
             'deeppink2',
             'firebrick4',
             'gold3',
             'slategray4',
             'magenta3',
             'darkorange3',
             'navajowhite3')
