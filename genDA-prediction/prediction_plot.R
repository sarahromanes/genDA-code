library(tidyverse)
load("prediction.RData")

## There were MANY plots produced from this data, all you need to do is change the filtering based on covairance, n, m, and sep to get the plot you want ##

df <- results.prediction %>% filter(covariance == "Unequal" & n == 100 & m == 50 & sep!= 0.75)

ggplot(df, aes(x = method, y = value, fill = method)) + 
  geom_boxplot() + 
  facet_wrap(sep~ ., scales = "free_y", nrow = 2)  +
  theme(axis.title = element_text(size = 18), 
        axis.text = element_text(size = 18), 
        axis.text.x = element_text(size = 16, angle = 90, hjust = 1), 
        legend.position = "none", 
        strip.text = element_text(size = 18),
        plot.title = element_text(size = 20)) +
  scale_color_brewer(palette = "Dark2") + 
  ylab("Error") +
  xlab("Method")

