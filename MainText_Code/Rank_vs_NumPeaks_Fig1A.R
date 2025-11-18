data = read.csv("Npeaks_vs_Rank_1A.csv")    #Csv file included in MainText_Data folder

custom_colors <- c("gen" = "red", "epi" = "deepskyblue")
custom_shapes <- c("gen" = 21, "epi" = 23)

df <- data %>%
  group_by(npeaks, cond) %>%
  mutate(npeaks_jittered = npeaks + 0.05) %>%
  ungroup()

#triangle <- data.frame(x = c(1, 1, 5),y = c(1, 5, 5))

tiff("Figures/Rank_vs_npeak.tiff", width = 5, height = 5, units = "in", res = 300)
ggplot(df, aes(x = npeaks_jittered, y = mean, color = cond, shape = cond)) +
  #geom_point(size = 4,color="black") + geom_point(size=2)+
  #geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1,linewidth = 1.0) +
  #scale_color_manual(values = custom_colors) + 
  #scale_shape_manual(values = custom_shapes) +
  annotate("segment", x = 1, xend = 7, y = 1, yend = 4,color = "gray50", linetype = "dashed", linewidth = 1)+
  #geom_polygon(data = triangle, aes(x = x, y = y), fill = "lightpink",inherit.aes = FALSE) +
  annotate("segment", x = 1, xend = 5, y = 1, yend = 5,color = "gray50", linetype = "solid", linewidth = 1)+
  geom_errorbar(aes(ymin = lower, ymax = upper,color=cond), width = 0.1, linewidth=1) + 
  geom_point(size = 4, stroke = 1.5, color = "black",aes(fill=cond)) + 
  # Set colors and shapes
  scale_fill_manual(values = custom_colors) + 
  scale_shape_manual(values = custom_shapes) +
  scale_color_manual(values = custom_colors) +
  theme_classic()+theme(text = element_text(size = 26),legend.position="none")+
  #theme_classic()+theme(text = element_text(size = 26))+
  xlab("Number of Peaks") + ylab("Mean Rank")+
  scale_x_continuous(breaks = 1:7)+
  scale_y_reverse()+theme(legend.title = element_blank(),legend.position.inside = c(1, 1),legend.justification = c(1, 0.5))
dev.off()
