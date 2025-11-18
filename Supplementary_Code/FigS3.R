#GlobalPeak Probability vs Npeaks
data = read.csv("Npeaks_vs_ProbGlobal_S3.csv")

custom_colors <- c("gen" = "red", "epi" = "deepskyblue")
custom_shapes <- c("gen" = 21, "epi" = 23)
line_data <- data.frame(x = seq(1, 7, length.out = 500),y = 1 /(seq(1, 7, length.out = 500)))

df <- data %>%
  group_by(npeaks, cond) %>%
  mutate(npeaks_jittered = if_else(cond == "gen", npeaks + 0.035, npeaks)) %>%
  ungroup()

tiff("Supplementary/Fig1/GlobalPeak_vs_npeak.tiff", width = 6, height = 6, units = "in", res = 300)
ggplot(df, aes(x = npeaks_jittered, y = mean, color = cond, shape = cond)) +
  geom_line(data = line_data, aes(x = x, y = y), inherit.aes = FALSE,linetype = "dashed", color = "gray50", linewidth = 1)+
  geom_errorbar(aes(ymin = lower, ymax = upper,color=cond), width = 0.1, linewidth=1) + 
  geom_point(size = 4, stroke = 1.5, color = "black",aes(fill=cond)) + 
  scale_fill_manual(values = custom_colors) + 
  scale_shape_manual(values = custom_shapes) +
  scale_color_manual(values = custom_colors) +
  theme_classic()+theme(text = element_text(size = 26),legend.position="none")+
  xlab("Number of Peaks") + ylab("Probability Global Peak")+
  scale_x_continuous(breaks = 1:7)
dev.off()

#MeanFitness vs Npeaks
data = read.csv("Npeaks_vs_ProbGlobal_S3.csv")

custom_colors <- c("gen" = "red", "epi" = "deepskyblue")
custom_shapes <- c("gen" = 21, "epi" = 23)

df <- data 

tiff("Supplementary/Fig1/Fitness_vs_npeak.tiff", width = 6, height = 6, units = "in", res = 300)
ggplot(df, aes(x = npeaks, y = meanfit, color = cond, shape = cond)) +
  geom_errorbar(aes(ymin = lower_fit, ymax = upper_fit,color=cond), width = 0.1, linewidth=1) + 
  geom_point(size = 4, stroke = 1.5, color = "black",aes(fill=cond)) + 
  scale_fill_manual(values = custom_colors) + 
  scale_shape_manual(values = custom_shapes) +
  scale_color_manual(values = custom_colors) +
  theme_classic()+theme(text = element_text(size = 26),legend.position="none")+
  xlab("Number of Peaks") + ylab("Mean Fitness")+
  scale_x_continuous(breaks = 1:7)
dev.off()
