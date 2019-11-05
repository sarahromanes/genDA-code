
library(tidyverse)

                       #----------------- Procrutes Plots ---------------#

load("Bernoulli_full.RData")
load("NB_full.RData")
load("mixed_genDA.RData")


# Bernoulli #

df <- results.bern %>% filter(parameter == "mU" | parameter == "mLambda")

ggplot(df, aes(x = method, y = value, fill = method)) +
  geom_boxplot() +
  facet_wrap(parameter ~ m +n, scales = "free_y", nrow = 2) +
  theme_bw() + 
  ylab("Procrustes error") + 
  theme(text = element_text(size=18), 
        legend.position = "none",  
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())

# Negative Binomial #

df <- results.NB %>% filter(parameter == "mU" | parameter == "mLambda")

ggplot(df, aes(x = method, y = value, fill = method)) +
  geom_boxplot() +
  facet_wrap(parameter ~ m +n, scales = "free_y", nrow = 2) +
  theme_bw() + 
  ylab("Procrustes error") + 
  theme(text = element_text(size=18), 
        legend.position = "none",  
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())


# Mixed Results #

df <- results.mixed.genDA %>% filter(parameter == "mU" | parameter == "mLambda")

ggplot(df, aes(x = method, y = value, fill = method)) +
  geom_boxplot() +
  facet_wrap(parameter ~ m +n, scales = "free_y", nrow = 2) +
  theme_bw() + 
  ylab("Procrustes error") + 
  theme(text = element_text(size=18), 
        legend.position = "none",  
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())




                          #----------------- Time Plots ---------------#

df = results.bern %>% filter(parameter =="time" & n == 50) %>% group_by(n, m, method) %>% summarise(avg = mean(value))
ggplot(df, aes(x = m, y = avg, colour = method)) + 
  geom_line(size = 1.5, alpha = 0.6) + 
  geom_point(size = 4) + 
  theme_bw() +
  theme(axis.title = element_text(size = 24), 
        axis.text = element_text(size = 24), 
        axis.text.x = element_text(size = 24), 
        legend.text = element_text(size = 24), 
        legend.title = element_text(size = 24), 
        strip.text = element_text(size = 24),
        plot.title = element_text(size = 30),
        legend.position = "bottom") +
  ylab("Time (s)") +
  ggtitle("Binary GLLVM")




df = results.NB %>% filter(parameter =="time" & n == 50) %>% group_by(n, m, method) %>% summarise(avg = mean(value))
ggplot(df, aes(x = m, y = avg, colour = method)) + 
  geom_line(size = 1.5, alpha = 0.6) + 
  geom_point(size = 4) + 
  theme_bw() +
  theme(axis.title = element_text(size = 24), 
        axis.text = element_text(size = 24), 
        axis.text.x = element_text(size = 24), 
        legend.text = element_text(size = 24), 
        legend.title = element_text(size = 24), 
        strip.text = element_text(size = 24),
        plot.title = element_text(size = 30),
        legend.position = "bottom") +
  ylab("Time (s)") +
  ggtitle("Over-dispersed Count GLLVM")



