data=read.csv("Csv files for figures/EpigeneticRates_HeatMap.csv")
allowed_values <- c(0.01,0.03,0.1,0.3,1)
filtered_data <- data %>%
  filter(RateOn %in% allowed_values, RateOff %in% allowed_values)

tiff("Supplementary/Fig 4/EpigeneticRates_FitnessHeatMap.tiff", width = 6, height = 6, units = "in", res = 300)
ggplot(filtered_data, aes(x = factor(RateOn), y = factor(RateOff), fill = mean_fit)) +
  geom_tile(color = "black") +
  #geom_text(aes(label = round(mean, 2)), color = "white", size = 5) +  # value labels
  scale_fill_gradient(low = "white", high = "gray20",name = "Mean Fitness") +
  #scale_fill_viridis_c(option = "plasma")+
  theme_minimal() +  theme(text = element_text(size = 22),legend.title = element_text(hjust =0.5))+
  xlab("Probability ON")+ylab("Probability OFF")+coord_fixed()
dev.off()
