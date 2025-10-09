#Histogram of number of peaks - Main 
xaxis = c(1,2,3,4,5,6,7)
nlands = c(36, 128, 243, 288, 189, 84, 27)
df <- data.frame(xaxis, nlands)
df$freq <- df$nlands / sum(df$nlands)

tiff("Hist_npeaks.tiff", width = 6, height = 6, units = "in", res = 300)
ggplot(df, aes(x = factor(xaxis), y = nlands)) +
  geom_col(fill = "lightgray", color="black",linewidth=1) +
  theme_classic()+theme(text = element_text(size = 26),legend.position="none")+
  xlab("Number of peaks") + ylab("Frequency")
dev.off()

#Histogram for number of peaks - House of Cards
xaxis = c(1,2,3,4,5,6,7,8,9)
nlands = c(0,8,74,196,288,259,123,41,9)
df <- data.frame(xaxis, nlands)
df$freq <- df$nlands / sum(df$nlands)

tiff("Hist_npeaks_HoC.tiff", width = 6, height = 6, units = "in", res = 300)
ggplot(df, aes(x = factor(xaxis), y = nlands)) +
  geom_col(fill = "lightgray", color="black",linewidth=1) +
  theme_classic()+theme(text = element_text(size = 26),legend.position="none")+
  xlab("Number of peaks") + ylab("Frequency")
dev.off()
