library(stringr)
library(MASS)
library(doParallel)
library(ggplot2)
library(flock)
library(foreach)

#generate a set of a large number of landscapes.
nsim = 1000  #number of landscapes
npeaks=4     #if we want to generate landscapes with fixed number of peaks, else not of interest

L = 5; nodes = 2^L                   #number of loci


combinations <- expand.grid(rep(list(c(0, 1)), L))          
genotypes <- apply(combinations, 1, paste0, collapse = "")  #generating list of all possible genotypes of length L

landscapes = matrix(0,nrow=nodes,ncol=nsim)   #matrix to store fitness values of all genotypes
numpeaks = numeric(nsim)                      #vector to store number of peaks of each landscape
listpeaks = list()                            #vector to store locus of each peak of each landscape
sdhoc = 0.25                                  #standard deviation of pairwise epistasis
rankpeaks = list()

#while(i<nsim) #use if generating landscapes of a given ruggedness or number of peaks
for(i in 1:nsim)  #loop for generating nsim number of landscapes
{
  print(i)
  fit = numeric(nodes); fit[1] = 1  #vector to store fitness of each genotype. defining WT fitness as 1.0
  peaks = numeric(nodes)            #vector to store if a given genotype is a peak
  #assign fitness to each genotype
  for(z in 2:nodes)
  {
    g = genotypes[z]
    nmut = str_count(g, "1")
    positions <- which(strsplit(g, "")[[1]] == "1")
    
    if(length(positions) > 1)
    {
      s = rnorm(1,0,sdhoc)
    } else 
    {
      s = 0
    }
    fit[z] = 1 + s
    if(fit[z] < 0)
    {
      fit[z] = 0
    }
  }
  
  #estimate number of peaks
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
    if(fit[z]>max(fit[ind]))
    {
      true=1
    }
    peaks[z] = true
  }
  indpeaks = which(peaks==1)
  
  #if(length(indpeaks)==npeaks)
  #{
  #  i=i+1
  #  landscapes[,i] = fit
  #  numpeaks[i] = length(indpeaks)
  #  listpeaks[[i]] = indpeaks
  #}
  landscapes[,i] = fit
  numpeaks[i] = length(indpeaks)
  listpeaks[[i]] = indpeaks
  rankpeaks[[i]] = rank(-fit[indpeaks])
}
