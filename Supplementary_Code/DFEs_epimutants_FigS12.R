#Code to make a table of DFEs - 
#change the set of landscapes in the "landscapes" variable for different sets of landscapes.
#defined for a 5-locus landscape. parameters will need to be changed if that is no longer true! (specifically 16, and 32)
dfe = numeric(nsim*16)
count = 0
for(i in 1:nsim)
{
  fit = landscapes[,i]
  for(j in 1:32)
  {
    if(j%%2 == 1)
    {
      count = count+1
      dfe[count] = fit[j+1] - fit[j]
    }
  }
}

write.csv(dfe, "dfe_epi_epist.csv", row.names = FALSE)

#making this histogram - 
#DFEs of epigenetics
#for HoC model
data = read.csv("House of Cards/dfe_epi_HoC.csv")
tiff("House of Cards/dfe_epi_Hoc.tiff", width = 5, height = 5, units = "in", res = 300)
ggplot(data, aes(x = x)) +
  #geom_histogram(aes(y = after_stat(density)), fill = "deepskyblue", color = "black", bins = 10) +
  geom_histogram(fill = "deepskyblue", color = "black", bins = 50) +
  ylab("Number of epimutations")+xlab("Fitness effect")+
  theme_classic()+theme(text = element_text(size = 26),legend.position="none")+xlim(-1.5,1.5)
dev.off()

#For main model
data = read.csv("dfe_epi_epist.csv")
tiff("dfe_epi_epist.tiff", width = 5, height = 5, units = "in", res = 300)
ggplot(data, aes(x = x)) +
  #geom_histogram(aes(y = after_stat(density)), fill = "deepskyblue", color = "black", bins = 10) +
  geom_histogram(fill = "deepskyblue", color = "black", bins = 50) +
  ylab("Number of epimutations")+xlab("Fitness effect")+
  theme_classic()+theme(text = element_text(size = 26),legend.position="none")+xlim(-1.5,1.5)
dev.off()
