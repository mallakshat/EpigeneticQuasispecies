##Fig S13A - data generated using the main WF model
data = read.csv("Npeaks_vs_Rank_HoC_S13A.csv")
custom_colors <- c("gen" = "red", "epi" = "deepskyblue")
custom_shapes <- c("gen" = 21, "epi" = 23)
df <- data 

tiff("House of Cards/RankHoC_vs_npeaks_new.tiff", width = 5, height = 5, units = "in", res = 300)
ggplot(df, aes(x = npeaks, y = mean, color = cond, shape = cond)) +
  annotate("segment", x = 1, xend = 8, y = 1, yend = 4.5,color = "gray50", linetype = "dashed", linewidth = 1)+
  annotate("segment", x = 1, xend = 5, y = 1, yend = 5,color = "gray50", linetype = "solid", linewidth = 1)+
  geom_errorbar(aes(ymin = lower, ymax = upper,color=cond), width = 0.1, linewidth=1) + 
  geom_point(size = 4, stroke = 1.5, color = "black",aes(fill=cond)) + 
  scale_fill_manual(values = custom_colors) + 
  scale_shape_manual(values = custom_shapes) +
  scale_color_manual(values = custom_colors) +
  theme_classic()+theme(text = element_text(size = 26),legend.position="none")+
  xlab("Number of Peaks") + ylab("Mean Rank")+
  scale_x_continuous(breaks = 1:8, limits=c(1,8.2))+scale_y_reverse(limits=c(1,5),breaks=1:5)
dev.off()

#Fig S13B - I removed cases where we had only 1 or 2 points - data generated using the script for Fig 2 in MainText_Code.
data <- read.csv("RankImprovement_HoC_S13B.csv")
data <- data %>% filter(epinumpeaks5 != 6)
data <- data %>% filter(epinumpeaks5 != 2)

tiff("House of Cards/RankImprovement.tiff", width = 5, height = 5, units = "in", res = 300)
ggplot(data, aes(x = (5-epinumpeaks5), y = delrank, group=epinumpeaks5))+
  #geom_violin(trim = FALSE, fill = "steelblue", color = "black") +
  geom_jitter(width = 0.1, size = 2, alpha = 0.4) +  # Add jittered points
  #geom_boxplot(width = 0.1, fill = "white", color = "darkblue", outlier.shape = NA) +
  theme_classic() + theme(text = element_text(size = 24))+
  xlab("Peaks removed") + ylab("Rank improvement")+
  geom_hline(yintercept = 0, linetype = "dashed",color="black",linewidth = 0.5)+
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, color = "black", linewidth = 0.7)
dev.off()

#Fig S13C - data generated using script for Fig S6 in Supplementary_Code
data = read.csv("Heterogeneity_EffectiveLands_HoC_S13C.csv")

shapes <- c("NA"= 1,"deepskyblue"= 1,"red" = 1)
tiff("Effective landscapes/Correlation_Eff_Epi_HoC.tiff", width = 5, height = 5, units = "in", res = 300)
ggplot(data, aes(x = rank_gen, y = rank_epi)) +
  #geom_point(data = subset(df, is.na(fill_color)), shape = 1, size = 3,stroke=1) +  # Unfilled circles
  #geom_point(data = subset(df, !is.na(fill_color)), aes(fill = fill_color), shape = 21, color = "black", size = 3, stroke=1) +
  geom_smooth(method = "lm", color = "black", aes(group = 1), se=FALSE)+
  geom_point(color="black",size=3, shape=1)+
  #scale_shape_manual(values = shapes)+
  #scale_color_manual(values = c("gray","deepskyblue","red"))+
  #geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black",linewidth=0.75) + 
  theme_classic()+theme(text = element_text(size = 26),legend.position="none")+
  ylab("Rank (G+Epi)") + xlab("Rank (Effective G)")+
  scale_x_reverse() + scale_y_reverse()
dev.off()

summary(lm(data$rank_gen ~ data$rank_epi))
