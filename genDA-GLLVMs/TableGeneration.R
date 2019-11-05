## Tables 4.2, 4.3, and 4.4 ##

# You can filter by n and m to make reading easier as well. I know you can probably automate table generation to LaTeX in R but I just wrote it in manually into my .tex file... 

library(tidyverse)

load("Bernoulli_full.RData")
load("NB_full.RData")
load("mixed_genDA.RData")

results.bern %>% filter(parameter!= "mU", parameter!="mLambda", parameter != "time") %>% group_by(parameter,n, m, method, metric) %>% summarise(avg=mean(value)) 

results.NB %>% filter(parameter!= "mU", parameter!="mLambda", parameter != "time") %>% group_by(parameter,n, m, method, metric) %>% summarise(avg=mean(value)) 

results.mixed.genDA %>% filter(parameter!= "mU", parameter!="mLambda", parameter != "time") %>% group_by(parameter,n, m, method, metric) %>% summarise(avg=mean(value)) 
