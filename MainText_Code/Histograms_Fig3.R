##Histogram for genetics
xaxis = rev(c("01100","10110","11101","11010","01111"))
freq = c(5,159,4,105,227)/reps
fitval = c(1.94 ,1.90, 1.81, 1.70, 1.57)
fitvals = c(1,2,3,4,5)

df <- data.frame(xaxis, freq, fitvals, fitval)
df$xaxis <- factor(df$xaxis, levels = xaxis)
custom_labels <- paste0(df$xaxis, "\n", df$fitvals, "\n", df$fitval)

tiff("Figures/Freq_genetics.tiff", width = 6, height = 6, units = "in", res = 300)
ggplot(df, aes(x = xaxis, y = freq)) +
  geom_col(fill = "red", color = "black",linewidth = 1.5) +  # Black border for bars
  theme_classic() + xlab("Genotype, Rank and Fitness") + ylab("Frequency")+
  theme(text = element_text(size = 28),legend.position="none")+
  scale_x_discrete(labels = custom_labels) +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.1))
dev.off()  



#Histogram for Epigenetics
xaxis = rev(c("01100","10110","11101","11010","01111"))
freq = c(90,277,43,90,0)/reps
fitval = c(1.93 ,1.90, 1.81, 1.70, 1.56)

df <- data.frame(xaxis, freq, fitvals, fitval)
df$xaxis <- factor(df$xaxis, levels = xaxis)
custom_labels <- paste0(df$xaxis, "\n", df$fitvals, "\n", df$fitval)


tiff("Figures/Freq_epigenetics.tiff", width = 6, height = 6, units = "in", res = 300)
ggplot(df, aes(x = xaxis, y = freq)) +
  geom_col(fill = "deepskyblue", color = "black",linewidth = 1.5) +  # Black border for bars
  theme_classic() + xlab("Genotype, Rank and Fitness") + ylab("Frequency")+
  theme(text = element_text(size = 28),legend.position="none")+
  scale_x_discrete(labels = custom_labels)+
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.1))
dev.off()
