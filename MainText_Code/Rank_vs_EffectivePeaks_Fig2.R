#Effective fitness with epigenetics
solve_quadratic_positive <- function(a, b, c) {
  discriminant <- b^2 - 4*a*c
  
  if (discriminant < 0) {
    return("No real solutions")
  } else {
    x1 <- (-b + sqrt(discriminant)) / (2*a)
    x2 <- (-b - sqrt(discriminant)) / (2*a)
    
    positive_solutions <- c(x1, x2)[c(x1, x2) >= 0 & c(x1, x2) <= 1]
    
    if (length(positive_solutions) == 0) {
      return("No viable solutions")
    } else {
      return(max(positive_solutions))
    }
  }
}


epilandscapes = matrix(0,nrow=nodes,ncol=nsim)
epinumpeaks = numeric(nsim)
epilistpeaks = list()
zerocount = 0

peak_rem = list()                                    #list of rank of peak removed
for (i in 1:nsim) 
{
  print(i)
  fit = landscapes[,i]
  
  indepi = matrix(0, ncol = length(loci), nrow = nodes)  #vector to store which genotype is the same as an epigenetic switch from every genotype
  for(bb in 1:length(loci))
  {
    locnow = loci[bb]
    changes = numeric(nodes)
    for(z in 1:nodes)
    {
      g = genotypes[z]
      split <- unlist(strsplit(g, ""))
      split[loci[bb]] <- "1"
      join <- paste(split, collapse = "")
      index = which(genotypes==join)
      changes[z] = index
    }
    stopImplicitCluster()
    indepi[,bb] = changes
  }
  fitepi = fit[indepi]             #fitness of each epigenetic "on" state of each genotype
  
  for (j in 1:nodes) 
  {
    if(indepi[j]==j)
    {
      epilandscapes[j,i] = fit[j]
      next
    }
    g = genotypes[j]
    genfitness = fit[j]
    epifitness = fitepi[j]
    s = (fitepi[j]/fit[j])-1
    coeff_a = s
    coeff_b = (1-rateup) - ((1+s)*(1+ratedown))
    coeff_c = (1+s)*ratedown
    pequil = solve_quadratic_positive(coeff_a, coeff_b, coeff_c)
    #if(length(pequil)!=1)
    #{
    #  print(j)
    #  print(pequil)
    #}
    #print(pequil)
    if(pequil==0)
    {
      zerocount=zerocount+1
    }
    epilandscapes[j,i] = genfitness*pequil + epifitness*(1-pequil)
  }
  
  fit_eff = epilandscapes[,i]
  #estimate number of effective peaks
  peaks=numeric(nodes)
  for(z in 1:nodes)
  {
    true=0
    g = genotypes[z]
    neighbors = character()
    for (j in 1:L) 
    {
      neighbors <- c(neighbors, paste0(substr(g, 1, j-1), as.character(1-as.integer(substr(g, j, j))), substr(g, j+1, nchar(g))))
    }
    ind = which(genotypes %in% neighbors)
    if(fit_eff[z]>max(fit_eff[ind]))
    {
      true=1
    }
    peaks[z] = true
  }
  indpeaks = which(peaks==1)
  epinumpeaks[i] = length(indpeaks)
  epilistpeaks[[i]] = indpeaks
  
  miss = setdiff(listpeaks[[i]] , epilistpeaks[[i]])
  ind_rem <- if(length(miss) > 0) which(listpeaks[[i]] %in% miss) else integer(0)
  peak_rem[[i]] = rankpeaks[[i]][ind_rem]
}  

#Making plot
peaks_interest = 5
indices5 = which(numpeaks==peaks_interest)
nsim2 = length(indices5)


epinumpeaks5 = epinumpeaks[indices5]
delrank = rankpeak[,1] - rankpeak[,2]

df <- data.frame(delrank,epinumpeaks5)
# Calculate means for each epinumpeaks value
means <- df %>%
  group_by(epinumpeaks5) %>%
  summarise(mean_value = mean(delrank))

##I saved df as a csv file, included in MainText_Data. I used the following code to make the figure - 
data <- read.csv("Rank_Improvement_Fig2.csv")

ggplot(data, aes(x = (5-epinumpeaks5), y = delrank, group=epinumpeaks5))+
  #geom_violin(trim = FALSE, fill = "steelblue", color = "black") +
  geom_jitter(width = 0.1, size = 2, alpha = 0.4) +  # Add jittered points
  #geom_boxplot(width = 0.1, fill = "white", color = "darkblue", outlier.shape = NA) +
  theme_classic() + theme(text = element_text(size = 24))+
  xlab("Peaks removed") + ylab("Rank improvement")+
  geom_hline(yintercept = 0, linetype = "dashed",color="black",linewidth = 0.5)+
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, color = "black", linewidth = 0.7)
#stat_summary(fun = mean, geom = "point", shape = 16, size = 3, color = "red")
#stat_summary(fun = mean, geom = "text", vjust = -7, size = 8, aes(label = sprintf("%.2f", after_stat(y))))+
#geom_smooth(data = means, aes(x = (5-epinumpeaks5), y = mean_value), method = "lm",color = "black", group = 1, linewidth = 1)
dev.off()


