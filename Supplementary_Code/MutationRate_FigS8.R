#Mu vs GlobalPeakProbability
data=read.csv("Supplementary/Fig 4/Mu_vs_ProbGlobal.csv")

custom_colors <- c("gen" = "red", "epi" = "deepskyblue")
custom_shapes <- c("gen" = 21, "epi" = 23)

df <- data 

tiff("Supplementary/Fig 4/Mu_vs_GlobalPeak.tiff", width = 6, height = 6, units = "in", res = 300)
ggplot(df, aes(x = log10(mu), y = mean, color = cond, shape = cond)) +
  geom_vline(xintercept = -6,linetype="dashed",color="gray50",linewidth=1)+
  geom_errorbar(aes(ymin = lower, ymax = upper,color=cond), width = 0.05, linewidth=1) + 
  geom_point(size = 4, stroke = 1.5, color = "black",aes(fill=cond)) + 
  scale_fill_manual(values = custom_colors) + scale_shape_manual(values = custom_shapes) +
  scale_color_manual(values = custom_colors) +
  theme_classic()+theme(text = element_text(size = 26),legend.position="none")+
  xlab("Mutation Rate µ") + ylab("Probability Global Peak")+
  scale_x_continuous(breaks = c(-7,-6,-5),labels = c(expression(10^-7),expression(10^-6),expression(10^-5)))
dev.off()

#Mu vs Mean fitness
data=read.csv("Supplementary/Fig 4/Mu_vs_ProbGlobal.csv")

custom_colors <- c("gen" = "red", "epi" = "deepskyblue")
custom_shapes <- c("gen" = 21, "epi" = 23)

df <- data 

tiff("Supplementary/Fig 4/Mu_vs_fitness.tiff", width = 6, height = 6, units = "in", res = 300)
ggplot(df, aes(x = log10(mu), y = meanfit, color = cond, shape = cond)) +
  geom_vline(xintercept = -6,linetype="dashed",color="gray50",linewidth=1)+
  geom_errorbar(aes(ymin = lower_fit, ymax = upper_fit,color=cond), width = 0.05, linewidth=1) + 
  geom_point(size = 4, stroke = 1.5, color = "black",aes(fill=cond)) + 
  scale_fill_manual(values = custom_colors) + scale_shape_manual(values = custom_shapes) +
  scale_color_manual(values = custom_colors) +
  theme_classic()+theme(text = element_text(size = 26),legend.position="none")+
  xlab("Mutation Rate µ") + ylab("Mean fitness")+
  scale_x_continuous(breaks = c(-7,-6,-5),labels = c(expression(10^-7),expression(10^-6),expression(10^-5)))
dev.off()
