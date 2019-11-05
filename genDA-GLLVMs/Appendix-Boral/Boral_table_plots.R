library(tidyverse)

load("Bernoulli_boral.RData")
load("Mixed_boral.RData")


load("Bernoulli_full.RData")
load("NB_full.RData")
load("mixed_genDA.RData")

                                  ###------ Table of Results -------##

results.boral.bern %>% filter(parameter!= "mU", parameter!="mLambda", parameter != "time") %>% group_by(parameter,n, m, method, metric) %>% summarise(avg=mean(value)) 

results.boral.mixed %>% filter(parameter!= "mU", parameter!="mLambda", parameter != "time") %>% group_by(parameter,n, m, method, metric) %>% summarise(avg=mean(value)) 


## Get subset of genDA simulations that correspond to boral - only use first 20% of data as 100 trials for boral vs 500 trials for genDA 

df =results.mixed.genDA %>% group_by(parameter,n, m, method, metric) %>% filter(metric != "procrustes" & metric != "s" & method =="genDA") %>% filter((n==50& m==100) | (n==100 & m ==50))
df = df[1:2400,]  
df %>% summarise(avg = mean(value)) 

df =results.bern %>% group_by(parameter,n, m, method, metric) %>% filter(metric != "procrustes" & metric != "s" & method =="genDA") %>% filter((n==50& m==100) | (n==100 & m ==50))
df = df[1:1600,]  
df %>% summarise(avg = mean(value)) 

                                  #### -------  Timings ------------###

# Bernoulli - genDA vs gllvm vs boral #

times <-  results.bern %>% filter(parameter =="time") %>% filter((n==50& m==100) | (n==100 & m ==50)) 
times <- times[1:400,]

df_time = rbind(results.boral.bern %>% filter(parameter =="time"), times)

ggplot(df_time, aes(y = value, x = method, fill = method)) +
  geom_boxplot() + 
  facet_wrap(parameter ~m + n) +
  theme_bw()+
  scale_y_continuous(trans='log10') +
  theme(text = element_text(size = 18),
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 18)) +
  ylab("time (s)") +
  scale_fill_brewer(palette = "Dark2")


# Mixed - genDA vs boral #

times <- results.mixed.genDA %>% filter(parameter =="time") %>% filter((n==50& m==100) | (n==100 & m ==50)) 
times <- times[1:200,]

df_time = rbind(results.boral.bern %>% filter(parameter =="time"), times)

ggplot(df_time, aes(y = value, x = method, fill = method)) +
  geom_boxplot() + 
  facet_wrap(parameter ~m + n) +
  theme_bw()+
  scale_y_continuous(trans='log10') +
  theme(text = element_text(size = 18),
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 18)) +
  ylab("time (s)") +
  scale_fill_brewer(palette = "Dark2")

                                  ### -------- Procrustes Plots ------- ###

# Bernoulli - genDA vs gllvm vs boral #

df <- results.bern %>% filter(parameter == "mU" | parameter =="mLambda")%>% filter((n==50& m==100) | (n==100 & m ==50))
df <- df[1:800,]

procrustes <- rbind(results.boral %>% filter(parameter == "mU"| parameter == "mLambda"), df)
ggplot(procrustes, aes(x = method, y = value, fill = method)) +
  geom_boxplot() +
  facet_wrap(parameter ~ m +n, scales = "free_y", nrow = 2) +
  theme_bw() + 
  ylab("Procrustes error") + 
  theme(text = element_text(size=18), 
        legend.position = "none",  
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 18)) +
  scale_fill_brewer(palette = "Dark2")

# Mixed - genDA vs boral #

df <- results.mixed.genDA %>% filter(parameter == "mU" | parameter =="mLambda")%>% filter((n==50& m==100) | (n==100 & m ==50))
df <- df[1:400,]

procrustes <- rbind(results.boral.mixed %>% filter(parameter == "mU"| parameter == "mLambda"), df)
ggplot(procrustes, aes(x = method, y = value, fill = method)) +
  geom_boxplot() +
  facet_wrap(parameter ~ m +n, scales = "free_y", nrow = 2) +
  theme_bw() + 
  ylab("Procrustes error") + 
  theme(text = element_text(size=18), 
        legend.position = "none",  
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 18)) +
  scale_fill_brewer(palette = "Dark2")
