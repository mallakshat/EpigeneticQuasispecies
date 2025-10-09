##The same main code as WFmodel.R was used to generate the csv file used here, with mutation rates (u) or switching rates on and off changed. The following code was used to make the figures.
#Fig 4A
data=read.csv("Csv files for figures/Mu_vs_Rank.csv")

custom_colors <- c("gen" = "red", "epi" = "deepskyblue")
custom_shapes <- c("gen" = 21, "epi" = 23)

df <- data %>%
  group_by(mu, cond) %>%
  #mutate(npeaks_jittered = npeaks + runif(1, -0.1, 0.1)) %>%
  ungroup()

tiff("Mu_vs_rank.tiff", width = 6, height = 6, units = "in", res = 300)
ggplot(df, aes(x = log10(mu), y = mean, color = cond, shape = cond)) +
  #geom_point(size = 4,color="black") + geom_point(size=2)+
  #geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1,linewidth = 1.0) +
  #scale_color_manual(values = custom_colors) + 
  #scale_shape_manual(values = custom_shapes) +
  geom_vline(xintercept = -6,linetype="dashed",color="gray50",linewidth=1)+
  geom_errorbar(aes(ymin = lower, ymax = upper,color=cond), width = 0.05, linewidth=1) + 
  geom_point(size = 4, stroke = 1.5, color = "black",aes(fill=cond)) + 
  # Set colors and shapes
  scale_fill_manual(values = custom_colors) + 
  scale_shape_manual(values = custom_shapes) +
  scale_color_manual(values = custom_colors) +
  theme_classic()+theme(text = element_text(size = 26),legend.position="none")+
  xlab("Mutation Rate µ") + ylab("Mean Rank")+
  scale_x_continuous(breaks = c(-7,-6,-5),labels = c(expression(10^-7),expression(10^-6),expression(10^-5)))+
  scale_y_reverse()
dev.off()

#Fig 4B
data <- read.csv("Csv files for figures/EpigeneticRates_Diagonal.csv")

tiff("Figures/EpigeneticRates_Rank_new.tiff", width = 5, height = 5, units = "in", res = 300)
ggplot(data, aes(x = log10(RateOn), y = mean)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1,linewidth=1) +
  geom_point(size=4,color="black") + geom_point(size=2.5,color="white")+
  geom_hline(yintercept = 3.17, linetype = "dashed", color = "red",linewidth=1.25)+ 
  theme_classic() + xlab("Rate On")+ylab("Mean Rank")+
  scale_x_continuous(breaks = c(0,-1,-2,-3,-4,-5),labels = c(expression(1),expression(10^-1),expression(10^-2),expression(10^-3),expression(10^-4),expression(10^-5)))+
  theme(text = element_text(size = 26))+scale_y_reverse()
dev.off()

#Fig 4C
data=read.csv("Csv files for figures/EpigeneticRates_HeatMap.csv")
allowed_values <- c(0.01,0.03,0.1,0.3,1)
filtered_data <- data %>%
  filter(RateOn %in% allowed_values, RateOff %in% allowed_values)

tiff("Figures/EpigeneticRatesHeatMap.tiff", width = 6, height = 6, units = "in", res = 300)
ggplot(filtered_data, aes(x = factor(RateOn), y = factor(RateOff), fill = mean)) +
  geom_tile(color = "black", linewidth=0.7) +
  #geom_text(aes(label = round(mean, 2)), color = "white", size = 5) +  # value labels
  scale_fill_gradient(low = "gray20", high = "white",name = "Mean Rank",breaks = scales::pretty_breaks(n = 3)) +
  #scale_fill_viridis_c(option = "plasma")+
  theme_minimal() +  theme(text = element_text(size = 22),legend.title = element_text(hjust =0.5))+
  xlab("Probability ON")+ylab("Probability OFF")+coord_fixed()+
  guides(fill = guide_colourbar(reverse = TRUE))
dev.off()

