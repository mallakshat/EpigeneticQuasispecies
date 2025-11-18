data <- read.csv("EpigeneticRates_Diagonal_S10.csv")

tiff("Supplementary/Fig 4/EpigeneticRates_Fitness.tiff", width = 6, height = 6, units = "in", res = 300)
ggplot(data, aes(x = log10(RateOn), y = mean_fit)) +
  geom_errorbar(aes(ymin = lower_fit, ymax = upper_fit), width = 0.1,linewidth=1) +
  geom_point(size=4,color="black") + geom_point(size=2.5,color="white")+
  geom_hline(yintercept = 1.89, linetype = "dashed", color = "red",linewidth=1.25) +
  theme_classic() + xlab("Rate On")+ylab("Mean Fitness")+
  scale_x_continuous(breaks = c(0,-1,-2,-3,-4,-5),labels = c(expression(1),expression(10^-1),expression(10^-2),expression(10^-3),expression(10^-4),expression(10^-5)))+
  theme(text = element_text(size = 26))
dev.off()
