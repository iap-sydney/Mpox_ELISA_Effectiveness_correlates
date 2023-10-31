#Driver.R Runs all scripts required to recreate analysis in the article 

#Run Setup script
source('Setup.R')

# Process raw data into usable format
source('Processing/01_process-eff-data.R')
source('Processing/02_process-immune-data.R')

#Analyse Data
Analysis=TRUE

if (Analysis){
  source('Analysis/01_Efficacy-fit.R')
  source('Analysis/02_immunogenicity-fit-2.R')
  source('Analysis/03_logistic-curve-fit.R')
  source('Analysis/04_Decay-model-fit.R')
  source('Analysis/05_Efficacy-decay.R')
}

#Generate Figures
Figures=TRUE

if (Figures){
  source('Analysis/06_Figure-3.R')
}