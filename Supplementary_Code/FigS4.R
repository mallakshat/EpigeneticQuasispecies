#Change for 4 or 6 peaked case. To generate the dataset, please refer to the script for heterogeneity for Fig 1B in MainText_Code
data = read.csv("Heterogeneity_4peaks_S4A.csv")
##For the 6peaked case, use this instead - 
#data = read.csv("Heterogeneity_6peaks_S4B.csv")

shapes <- c("darkgray"= 1,"deepskyblue"= 16,"red" = 16)
tiff("Figures/Heterogeneity_Rank.tiff", width = 5, height = 5, units = "in", res = 300)
ggplot(data, aes(x = rank_gen, y = rank_epi, shape = col)) +
  #geom_point(data = subset(df, is.na(fill_color)), shape = 1, size = 3,stroke=1) +  # Unfilled circles
  #geom_point(data = subset(df, !is.na(fill_color)), aes(fill = fill_color), shape = 21, color = "black", size = 3, stroke=1) +
  geom_point(color="black",size=3)+
  scale_shape_manual(values = shapes)+
  #  scale_color_manual(values = c("gray","deepskyblue","red"))+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black",linewidth=0.75) + 
  theme_classic()+theme(text = element_text(size = 26),legend.position="none")+
  ylab("Rank (G+Epi)") + xlab("Rank (G only)")+
  scale_x_reverse() + scale_y_reverse()
dev.off()
