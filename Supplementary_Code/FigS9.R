data <- read.csv("Epigeneticrates_NewRatio_S9.csv")

tiff("Fig 4/EpigeneticDiagonal_NewRatio.tiff", width = 5, height = 5, units = "in", res = 300)
ggplot(data, aes(x = log10(RateOn), y = mean)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1,linewidth=1) +
  geom_point(size=4,color="black") + geom_point(size=2.5,color="white")+
  geom_hline(yintercept = 3.17, linetype = "dashed", color = "red",linewidth=1.25)+ 
  theme_classic() + xlab("Rate On")+ylab("Mean Rank")+
  scale_x_continuous(breaks = c(0,-1,-2,-3,-4,-5,-6),labels = c(expression(1),expression(10^-1),expression(10^-2),expression(10^-3),expression(10^-4),expression(10^-5),expression(10^-6)))+
  theme(text = element_text(size = 26))+scale_y_reverse()+ylim(3.5,2.5)
dev.off()
