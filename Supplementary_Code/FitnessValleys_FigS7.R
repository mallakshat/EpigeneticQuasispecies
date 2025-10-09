delta_gen = numeric(1)
delta_epi = numeric(1)
s_gen = numeric(1)
s_epi = numeric(1)
count = 0

for(xxx in 1:nsim2) #checking over all landscapes or x-peaked ones?
{
  i = indices5[xxx]
  fit = landscapes[,i]
  epifit = epilandscapes[,i]
  for(j in 1:numpeaks[[i]])
  {
    if(rankpeaks[[i]][j] == 1)
    {
      next
    }
    genotype_peak = genotypes[listpeaks[[i]][j]]
    
    chars <- strsplit(genotype_peak, "")[[1]]
    n <- length(chars)
    combs <- combn(n, 2)
    neighbors <- apply(combs, 2, function(pos) {
      new_chars <- chars
      new_chars[pos] <- ifelse(new_chars[pos] == "0", "1", "0")
      paste(new_chars, collapse = "")
    })
    neighbor_indices <- match(neighbors, genotypes)
    neighbor_indices
    
    for(k in 1:length(neighbor_indices))
    {
      if(fit[neighbor_indices[k]] > fit[listpeaks[[i]][j]]) #There are 2 possible valleys for every pair of genotypes?
      {
        c1 <- strsplit(genotype_peak, "")[[1]]
        c2 <- strsplit(neighbors[k], "")[[1]]
        diff_pos <- which(c1 != c2)
        intermediates <- sapply(diff_pos, function(pos) {
          tmp <- c1
          tmp[pos] <- c2[pos]
          paste(tmp, collapse = "")
        })
        intermediates_indices <- match(intermediates, genotypes)
        intermediates_indices
        
        count = count+1
        s_gen[count] = (fit[neighbor_indices[k]] - fit[listpeaks[[i]][j]])/fit[listpeaks[[i]][j]]
        delta_gen[count] = (fit[listpeaks[[i]][j]] - fit[intermediates_indices[1]])/fit[listpeaks[[i]][j]]
        s_epi[count] = (epifit[neighbor_indices[k]] - epifit[listpeaks[[i]][j]])/epifit[listpeaks[[i]][j]]
        delta_epi[count] = (epifit[listpeaks[[i]][j]] - epifit[intermediates_indices[1]])/epifit[listpeaks[[i]][j]]
        
        count = count+1
        s_gen[count] = (fit[neighbor_indices[k]] - fit[listpeaks[[i]][j]])/fit[listpeaks[[i]][j]]
        delta_gen[count] = (fit[listpeaks[[i]][j]] - fit[intermediates_indices[2]])/fit[listpeaks[[i]][j]]
        s_epi[count] = (epifit[neighbor_indices[k]] - epifit[listpeaks[[i]][j]])/epifit[listpeaks[[i]][j]]
        delta_epi[count] = (epifit[listpeaks[[i]][j]] - epifit[intermediates_indices[2]])/epifit[listpeaks[[i]][j]]
        
      }
    }
    
  }
}
delta_gen = -delta_gen
delta_epi = -delta_epi

df = data.frame(delta_gen, s_gen, delta_epi, s_epi)
write.csv(df, "Valleys_Delta_s_5peaked.csv", row.names = FALSE)

##Using this csv file, I make the figure following - 
data = read.csv("Valleys_Delta_s_5peaked.csv")

#First, deleterious steps size
df_long <- data %>%
  select(delta_gen, delta_epi) %>%
  pivot_longer(cols = everything(),names_to = "Variable",values_to = "Value")
df_long$Variable <- factor(df_long$Variable, levels = c("delta_gen", "delta_epi"))

summary_df <- df_long %>%
  group_by(Variable) %>%
  summarise(mean = mean(Value, na.rm = TRUE),se = sd(Value, na.rm = TRUE) / sqrt(n()),n = n()) %>%
  mutate(ci = 1.96 * se)

summary_df$group <- factor(summary_df$Variable,levels = c("delta_gen", "delta_epi"),labels = c("G only", "G+E"))

tiff("delta_valleys.tiff", width = 5, height = 5, units = "in", res = 300)
ggplot(summary_df, aes(x = group, y = -mean, fill = group)) +
  geom_col(width = 0.6, color = "black") +
  geom_errorbar(aes(ymin = -mean - ci, ymax = -mean + ci),width = 0.05, linewidth = 1) +
  scale_fill_manual(values = c("red", "deepskyblue")) +
  labs(x = "", y = "δ") +
  theme_classic(base_size = 26) +theme(legend.position = "none")
dev.off()

##Selective coefficients - s
df_long <- data %>%
  select(s_gen, s_epi) %>%
  pivot_longer(cols = everything(),names_to = "Variable",values_to = "Value")
df_long$Variable <- factor(df_long$Variable, levels = c("s_gen", "s_epi"))

summary_df <- df_long %>%
  group_by(Variable) %>%
  summarise(mean = mean(Value, na.rm = TRUE),se = sd(Value, na.rm = TRUE) / sqrt(n()),n = n()) %>%
  mutate(ci = 1.96 * se)
summary_df$group <- factor(summary_df$Variable,levels = c("s_gen", "s_epi"),labels = c("G only", "G+E"))

tiff("S_valleys.tiff", width = 5, height = 5, units = "in", res = 300)
ggplot(summary_df, aes(x = group, y = mean, fill = group)) +
  geom_col(width = 0.6, color = "black") +
  geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci),width = 0.05, linewidth = 1) +
  scale_fill_manual(values = c("red", "deepskyblue")) +
  labs(x = "", y = "s") +
  theme_classic(base_size = 26) +theme(legend.position = "none")
dev.off()

