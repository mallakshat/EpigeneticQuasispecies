#Fig S5A - epigenetics helps reduce the number of peaks
xaxis = (c(1,2,3,4,5))
freq = c(0,6,46,76,61)/189
df <- data.frame(xaxis, freq)

tiff("Supplementary/Fig 2/Freq_effectivePeaks.tiff", width = 5, height = 5, units = "in", res = 300)
ggplot(df, aes(x = xaxis, y = freq)) +
  geom_col(fill = "lightgreen", color = "black",linewidth = 1.5) +  # Black border for bars
  theme_classic() + xlab("Effective number of Peaks") + ylab("Frequency")+
  theme(text = element_text(size = 26),legend.position="none")+
  #scale_x_discrete(labels = custom_labels)+
  scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1))
dev.off()


#Fig S5B - rank of peak being removed.
data = read.csv("RankOfPeakRemoved.csv")
xvals <- c(1, 2, 3, 4, 5)
df <- data.frame(x = xvals, y = data$x/sum(data$x))

tiff("RankOfPeakRemoved.tiff", width = 5, height = 5, units = "in", res = 300)
ggplot(df, aes(x = factor(x), y = y)) +
  geom_col(fill = "lightgray", color="black",linewidth=1) +
  theme_classic()+theme(text = element_text(size = 26),legend.position="none")+
  xlab("Rank of peak removed") + ylab("Fraction")
dev.off()
