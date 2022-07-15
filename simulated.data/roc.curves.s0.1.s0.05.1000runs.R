##Thu  1 Aug 2019 10:40:54 BST
## R code used to assess power of 3P-CLR to detect selection in the context of chimpanzee demography (Aim 1). 
## load the pROC package (Robin et al. 2011)
library(pROC)
## directory containing the merged files (1000neutral + 1000s0.1), (1000neutral + 1000s0.05)
require(data.table)
## Assess power of 3P-CLR to detect selection as opposed to neutrality for the first selective coefficient modelled (s=0.1), using the Central/Eastern/Nigeria-Cameroon 3 population tree.

## 1a. load in the data for s=0.1
simulation_0.1_scores<-fread("merged.neutral.s0.1.likelihoodscores.1000runs.txt")
head(simulation_0.1_scores)
## add informative headers to the columns of the input file
colnames(simulation_0.1_scores)[1] <- "score_3pclr"
colnames(simulation_0.1_scores)[2] <- "run_identifier"

## Input file (merged.neutral.s0.1.likelihoodscores.1000runs.txt) contains the 3P-CLR scores from the central window of each simulation run in column 1, where each row corresponds to an independent simulation run, in total data from 1000 neutral and 1000 s_0.1 runs was used. 
## Column 2 of input file lists each simulation run's corresponding identifier - ie whether the 3P-CLR score listed corresponds to a neutral or selection simulation.   

head(simulation_0.1_scores)
## check variable type of the columns in the data table
str(simulation_0.1_scores)
## the predictor variable, ie the 3P-CLR scores (in column 1) needs to be a numeric continuous variable for input into pROC
## the response variable, ie the run identifier (in column 2) needs to be binary - as it can only have 2 levels, here this corresponds to either a neutral (encoded as 0) or selection (1) simulation


## estimate the ROC curve and calculate the area under the curve (AUC)
## values for AUC are used to compare ROC curves & hence compare accuracy of different models
roc(simulation_0.1_scores$run_identifier, simulation_0.1_scores$score_3pclr)

#Call:
 # roc.default(response = simulation_0.1_scores$run_identifier,     predictor = simulation_0.1_scores$score_3pclr)

#Data: simulation_0.1_scores$score_3pclr in 1000 controls (simulation_0.1_scores$run_identifier neutral) < 1000 cases (simulation_0.1_scores$run_identifier selection).
#Area under the curve: 0.9959

roc_s0.1 <- plot.roc(simulation_0.1_scores$run_identifier, simulation_0.1_scores$score_3pclr,main="Power comparison", col="darkblue")


##################

## 1b. load in the data for s=0.05
simulation_0.05_scores<-fread("merged.neutral.s0.05.likelihoodscores.1000runs.txt")
colnames(simulation_0.05_scores)[1] <- "score_3pclr"
colnames(simulation_0.05_scores)[2] <- "run_identifier"
head(simulation_0.05_scores)
roc(simulation_0.05_scores$run_identifier, simulation_0.05_scores$score_3pclr)

#Call:
 # roc.default(response = simulation_0.05_scores$run_identifier,     predictor = simulation_0.05_scores$score_3pclr)

#Data: simulation_0.05_scores$score_3pclr in 1000 controls (simulation_0.05_scores$run_identifier neutral) < 1000 cases (simulation_0.05_scores$run_identifier selection).
#Area under the curve: 0.9942

## 2. plot ROC curves for each selective coefficient

roc_s0.1 <- plot.roc(simulation_0.1_scores$run_identifier, simulation_0.1_scores$score_3pclr, col="darkblue")

## add 2nd plot to original plot
roc_s0.05<-lines.roc(simulation_0.05_scores$run_identifier, simulation_0.05_scores$score_3pclr, col="darkcyan")


## add legend to plot, corresponding to coloured curve for each s

legend("bottomright", legend=c("s=0.1", "s=0.05"), col=c("darkblue", "darkcyan"), lwd=3) 

##########################################################################################################

## 3. focus on top left-hand quadrant - visualise differences in power for each s at higher resolution

roc_s0.1_quadrant <- plot.roc(simulation_0.1_scores$run_identifier, simulation_0.1_scores$score_3pclr, col="darkblue", xlim=c(1.0,0.8), ylim=c(0.8,1.0))

## add 2nd plot to original plot
roc_s0.05_quadrant <- lines.roc(simulation_0.05_scores$run_identifier, simulation_0.05_scores$score_3pclr, col="darkcyan")
