### Data

- final_all_neutralruns_nc_withheader_all7elements.txt

-- input data to calculate F3 drift times for simulated data under neutrality


### To generate Figure S1: Power of 3P-CLR in the ancestral central-eastern population. Each ROC curve was generated from 1000 neutral simulations and 1000 selection simulations for each s (s=0.05, 0.1).

2 merged.neutral.s*.likelihoodscores.1000runs.txt files = 3P-CLR likelihood scores (3P-CLR calculated on simulated data) with identifier (as to whether score corresponds to a neutral simulation, or a simulation under positive selection)


- merged.neutral.s0.1.likelihoodscores.1000runs.txt

-- contains 3P-CLR scores from the central window of each simulation run in column 1, where each row corresponds to an independent simulation run, in total data from 1000 neutral and 1000 s_0.1 runs was used. 
Column 2 of input file lists each simulation run's corresponding identifier - ie whether the 3P-CLR score listed corresponds to a neutral or selection simulation.   
-- simulations generated using slim


- merged.neutral.s0.05.likelihoodscores.1000runs.txt

-- equivalent format but for selection coefficient s=0.05


- roc.curves.s0.1.s0.05.1000runs.R
-- code to generate ROC


### Code to generate Figure S11: Distribution of 3P-CLR scores for data simulated under neutrality (A) and positive selection with selection coefficients of 0.05 (B) and 0.1 (C), 1000 replicates were generated in each case.

```
sim_0.1<-fread("merged.neutral.s0.1.likelihoodscores.1000runs.txt")
sim_0.05<-fread("merged.neutral.s0.05.likelihoodscores.1000runs.txt")


neutral<-subset(sim_0.1, V2=="neutral")
n<-ggplot(neutral, aes(x=V1)) + geom_histogram(binwidth=.5) + labs(y="Frequency", x = "3P-CLR score") + ggtitle("Neutral simulations")
# distribution of 3PCLR scores for neutral simulated data

s_0.1<-subset(sim_0.1, V2=="selection")
x<-ggplot(s_0.1, aes(x=V1)) + geom_histogram(binwidth=.5) + labs(y="Frequency", x = "3P-CLR score") + ggtitle("Selection simulations (s=0.1)")
# add lines for the CIs 

s_0.05<-subset(sim_0.05, V2=="selection")
y<-ggplot(s_0.05, aes(x=V1)) + geom_histogram(binwidth=.5) + labs(y="Frequency", x = "3P-CLR score") + ggtitle("Selection simulations (s=0.05)")

# arrange on same page
library(ggpubr)
ggarrange(n, y, x, 
          labels = c("A", "B", "C"),
          ncol = 2, nrow = 2)
```
