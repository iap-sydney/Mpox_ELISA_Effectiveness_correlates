#Driver.R Runs all scripts required to recreate analysis in the article 

#Run Setup script
source('Setup.R')

#Create Output Folders
create_folders= !file.exists('Output/Figures')
if (create_folders){
  dir.create('Output/Data',recursive = TRUE)
  dir.create('Output/Samples',recursive = TRUE)
  dir.create('Output/Figure1',recursive = TRUE)
  dir.create('Output/Figure2',recursive = TRUE)
  dir.create('Output/Figure3',recursive = TRUE)
  dir.create('Output/Figure4',recursive = TRUE)
  dir.create('Output/Figures',recursive = TRUE)
}

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
  source('Analysis/06_Figure-1.R')
  source('Analysis/06_Figure-2.R')
  source('Analysis/06_Figure-3.R')
  source('Analysis/06_Figure-4.R')
  source('Analysis/07_Supplementary-Figures.R')
}
