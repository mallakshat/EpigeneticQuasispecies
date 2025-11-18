#define mutation rate
u = 1e-6
mut = u

#define epigenetic locus
loci = c(1)                                            #location which has an epigenetic switch
indepi = matrix(0, ncol = length(loci), nrow = nodes)  #vector to store which genotype is the same as an epigenetic switch from every genotype
for(bb in 1:length(loci))
{
  locnow = loci[bb]
  registerDoParallel(cores = 4)
  changes <- foreach(z = 1:nodes,.combine = c)%dopar%
    {
      g = genotypes[z]
      split <- unlist(strsplit(g, ""))
      split[loci[bb]] <- "1"
      join <- paste(split, collapse = "")
      index = which(genotypes==join)
      index
    }
  stopImplicitCluster()
  indepi[,bb] = changes
}

#Simulating adaptation
peaks_interest = 6
indices5 = which(numpeaks==peaks_interest)
nsim2 = length(indices5)

#data1
rd = c(0)
meanfit = matrix(0,nrow = nsim2, ncol = length(rd))
rankpeak = matrix(0,nrow = nsim2, ncol = length(rd))
entropy = matrix(0,nrow = nsim2, ncol = length(rd))
freqpeaks = matrix(0,nrow = nsim2, ncol = length(rd))
data1 = matrix(0,nrow = reps,ncol = nsim2)           #reps is initialized inside the loop, correct that.
for (count2 in (1:length(rd))) 
{
  #ratedown = rd[count2]                #rate of switching off of epigenetic state
  ratedown = rd[count2]
  rateup = ratedown*0.1
  for (count in 1:nsim2) 
  {
    print(count)
    ##WF model.
    k = 1e+5                              #carrying capacity
    fit = numeric(nodes)                  #variable to store fitness of each of the 2^L genotypes
    fit = landscapes[,indices5[count]]              #assigning fitness to fit from landscapes
    indpeaks=listpeaks[[indices5[count]]]           #storing indices of the fitness peaks
    indrankpeaks = rankpeaks[[indices5[count]]]
    fitepi = fit[indepi]                  #fitness of each epigenetic "on" state of each genotype
    
    reps = 200                            #number of simulations for one landscape
    registerDoParallel(cores = 48)
    #for(a in 1:reps)
    p <- foreach(a = 1:reps,.combine = 'rbind',.errorhandling = 'pass')%dopar%
      tryCatch({
        #p = numeric(reps)              #vector to store peak reached for every rep
        #initial conditions
        ini = numeric(nodes); ininum=0  #number of initial genotypes other than WT
        inis = sample(1:nodes,ininum)   #Note to self - change to sample from valleys, currently it is random
        ini[inis] = k/(ininum+1); ini[1] = k/(ininum+1)
        
        #initial epigenetic
        epi = numeric(nodes)            #number of switched individuals of each genotype
        
        time = 1e+4                     #number of generations
        n = matrix(0, nrow = nodes, ncol = time)  #matrix to store number of individuals of each genotype in every generation
        n[,1] = ini
        numepi = matrix(0, nrow = nodes, ncol = time)  #matrix to store number of individuals of each genotype in every generation
        
        for(i in 1:(time-1))
        {
          #epi[is.na(epi)] = 0
          #n[,i][is.na(n[,i])] = 0
          
          #Births
          fitbirth = fit*(n[,i]-epi) + fitepi*epi   #fitness for pbirth calculations based on impact of epigenetic switches
          pb = (fitbirth)/sum(fitbirth)       #probability of birth of each genotype
          nb = rbinom(nodes,k,pb)             #new births
          pbepi = (fitepi*epi)/(fitbirth)     #probability of descendant of a genotype being a mutant epiallele based of fitness differences
          pbepi[is.na(pbepi)] <- 0            #taking care of zero/zero scenarios
          epi = rbinom(nodes,nb,pbepi)        #new number of epigenetic mutants of each genotype after birth step
          
          #Epigenetic changes
          ons = rbinom(nodes,nb-epi,rateup)
          offs = rbinom(nodes,epi,ratedown)
          epi = epi+ons-offs
          
          #mutations
          nmut = rbinom(nodes,nb,mut)         #mutations in population
          indmut = which(nmut!=0)             #indices of genotypes chosen to mutate
          nmut[is.na(nmut)] <- 0 
          if(any(nmut!=0))
          {
            for(j in 1:length(indmut))
            {
              for(xyz in 1:nmut[indmut[j]])
              {
                x = indmut[j]             #which genotype is chosen for mutation
                y = genotypes[x]          #what's the exact genotype 
                #loc = which(strsplit(y, "")[[1]] == "0")
                loc = c(1:L)              #list of loci which can mutate 
                z <- unlist(strsplit(y, ""))
                posmut = sample(loc,1)    #position to mutate
                z[posmut] <- as.character(1 - as.numeric(z[posmut]))
                newmut <- paste(z, collapse = "")  #new genotype
                newind = which(genotypes==newmut)  #index of new genotype
                #does the individual chosen to mutate have an epigenetic on state
                probepi = epi[x]/nb[x] #likelihood of an individual from that genotype have an epigenetic change
                if(runif(1)<probepi)
                {
                  epi[x] = epi[x] - 1
                  epi[newind] = epi[newind] + 1
                }
                #changing genotype for chosen individual
                nb[newind] = nb[newind] + 1
                nb[x] = nb[x] - 1
                
              }
            }
          }
          
          #assigning to next generation
          n[,i+1] = nb
          lastepi=epi
          numepi[,i+1] = epi
          #check=0
          #win = which(nb==max(nb))
          #if(any(win==indpeaks))
          #{
          #  if(nb[win]>(0.999*k))
          #  {
          #    check=1
          #  }
          #}
          #if(check==1)
          #{
          #  break
          #}
        }
        
        #Outcome
        win = which(n[,time]==max(n[,time]))
        val = any(win==indpeaks)
        sorted_peaks = indpeaks[order(fit[indpeaks], decreasing = TRUE)]
        #peak_rank <- match(win, sorted_peaks)
        if(val==TRUE)
        {
          peak_rank = indrankpeaks[which(indpeaks==win)]
        } else
        {
          peak_rank = NA
        }
        
        
        finalfitness = fit*(n[,time]-epi) + fitepi*epi
        finalfit = sum(finalfitness)/sum(n[,time])
        c(peak_rank,finalfit)
      },error=function(e){
        return(c(NA, NA, paste("Error in rep", a, ":", e$message)))
      })
    stopImplicitCluster()
    
    meanfit[count,count2] = mean(p[,2])
  
    data1[,count] = p[,1]
    
    p = p[complete.cases(p), ]
    bstrapreps = (length(p[,1])) 
    bstraprank = numeric(bstrapreps)
    #bootstrap mean for reps number of times
    for(ii in 1:1000)
    {
      bstraprank[ii] = mean(sample(p[,1], bstrapreps, replace = TRUE))
    }
    rankpeak[count,count2] = mean(bstraprank)
    
    rankfreq = table(p[,1])/length(p[,1])
    entropy[count,count2] = -sum(rankfreq * log2(rankfreq))
    
    freqpeaks[count,count2] = length(p[,1])/reps
  }
  
}

#data2
rd = c(0.1)
meanfit = matrix(0,nrow = nsim2, ncol = length(rd))
rankpeak = matrix(0,nrow = nsim2, ncol = length(rd))
entropy = matrix(0,nrow = nsim2, ncol = length(rd))
freqpeaks = matrix(0,nrow = nsim2, ncol = length(rd))
data2 = matrix(0,nrow = reps,ncol = nsim2)
for (count2 in (1:length(rd))) 
{
  #ratedown = rd[count2]                #rate of switching off of epigenetic state
  ratedown = rd[count2]
  rateup = ratedown*0.1
  for (count in 1:nsim2) 
  {
    print(count)
    ##WF model.
    k = 1e+5                              #carrying capacity
    fit = numeric(nodes)                  #variable to store fitness of each of the 2^L genotypes
    fit = landscapes[,indices5[count]]              #assigning fitness to fit from landscapes
    indpeaks=listpeaks[[indices5[count]]]           #storing indices of the fitness peaks
    indrankpeaks = rankpeaks[[indices5[count]]]
    fitepi = fit[indepi]                  #fitness of each epigenetic "on" state of each genotype
    
    reps = 200                            #number of simulations for one landscape
    registerDoParallel(cores = 48)
    #for(a in 1:reps)
    p <- foreach(a = 1:reps,.combine = 'rbind',.errorhandling = 'pass')%dopar%
      tryCatch({
        #p = numeric(reps)              #vector to store peak reached for every rep
        #initial conditions
        ini = numeric(nodes); ininum=0  #number of initial genotypes other than WT
        inis = sample(1:nodes,ininum)   #Note to self - change to sample from valleys, currently it is random
        ini[inis] = k/(ininum+1); ini[1] = k/(ininum+1)
        
        #initial epigenetic
        epi = numeric(nodes)            #number of switched individuals of each genotype
        
        time = 1e+4                     #number of generations
        n = matrix(0, nrow = nodes, ncol = time)  #matrix to store number of individuals of each genotype in every generation
        n[,1] = ini
        numepi = matrix(0, nrow = nodes, ncol = time)  #matrix to store number of individuals of each genotype in every generation
        
        for(i in 1:(time-1))
        {
          #epi[is.na(epi)] = 0
          #n[,i][is.na(n[,i])] = 0
          
          #Births
          fitbirth = fit*(n[,i]-epi) + fitepi*epi   #fitness for pbirth calculations based on impact of epigenetic switches
          pb = (fitbirth)/sum(fitbirth)       #probability of birth of each genotype
          nb = rbinom(nodes,k,pb)             #new births
          pbepi = (fitepi*epi)/(fitbirth)     #probability of descendant of a genotype being a mutant epiallele based of fitness differences
          pbepi[is.na(pbepi)] <- 0            #taking care of zero/zero scenarios
          epi = rbinom(nodes,nb,pbepi)        #new number of epigenetic mutants of each genotype after birth step
          
          #Epigenetic changes
          ons = rbinom(nodes,nb-epi,rateup)
          offs = rbinom(nodes,epi,ratedown)
          epi = epi+ons-offs
          
          #mutations
          nmut = rbinom(nodes,nb,mut)         #mutations in population
          indmut = which(nmut!=0)             #indices of genotypes chosen to mutate
          nmut[is.na(nmut)] <- 0 
          if(any(nmut!=0))
          {
            for(j in 1:length(indmut))
            {
              for(xyz in 1:nmut[indmut[j]])
              {
                x = indmut[j]             #which genotype is chosen for mutation
                y = genotypes[x]          #what's the exact genotype 
                #loc = which(strsplit(y, "")[[1]] == "0")
                loc = c(1:L)              #list of loci which can mutate 
                z <- unlist(strsplit(y, ""))
                posmut = sample(loc,1)    #position to mutate
                z[posmut] <- as.character(1 - as.numeric(z[posmut]))
                newmut <- paste(z, collapse = "")  #new genotype
                newind = which(genotypes==newmut)  #index of new genotype
                #does the individual chosen to mutate have an epigenetic on state
                probepi = epi[x]/nb[x] #likelihood of an individual from that genotype have an epigenetic change
                if(runif(1)<probepi)
                {
                  epi[x] = epi[x] - 1
                  epi[newind] = epi[newind] + 1
                }
                #changing genotype for chosen individual
                nb[newind] = nb[newind] + 1
                nb[x] = nb[x] - 1
                
              }
            }
          }
          
          #assigning to next generation
          n[,i+1] = nb
          lastepi=epi
          numepi[,i+1] = epi
          #check=0
          #win = which(nb==max(nb))
          #if(any(win==indpeaks))
          #{
          #  if(nb[win]>(0.999*k))
          #  {
          #    check=1
          #  }
          #}
          #if(check==1)
          #{
          #  break
          #}
        }
        
        #Outcome
        win = which(n[,time]==max(n[,time]))
        val = any(win==indpeaks)
        sorted_peaks = indpeaks[order(fit[indpeaks], decreasing = TRUE)]
        #peak_rank <- match(win, sorted_peaks)
        if(val==TRUE)
        {
          peak_rank = indrankpeaks[which(indpeaks==win)]
        } else
        {
          peak_rank = NA
        }
        
        
        finalfitness = fit*(n[,time]-epi) + fitepi*epi
        finalfit = sum(finalfitness)/sum(n[,time])
        c(peak_rank,finalfit)
      },error=function(e){
        return(c(NA, NA, paste("Error in rep", a, ":", e$message)))
      })
    stopImplicitCluster()
    
    meanfit[count,count2] = mean(p[,2])
    
    data2[,count] = p[,1]
    
    p = p[complete.cases(p), ]
    bstrapreps = (length(p[,1])) 
    bstraprank = numeric(bstrapreps)
    #bootstrap mean for reps number of times
    for(ii in 1:1000)
    {
      bstraprank[ii] = mean(sample(p[,1], bstrapreps, replace = TRUE))
    }
    rankpeak[count,count2] = mean(bstraprank)
    
    rankfreq = table(p[,1])/length(p[,1])
    entropy[count,count2] = -sum(rankfreq * log2(rankfreq))
    
    freqpeaks[count,count2] = length(p[,1])/reps
    
  }
  
}


##Identifying statistically different landscapes.
rank_gen = numeric(nsim2)
rank_epi = numeric(nsim2)
col = character(nsim2)

for(i in 1:nsim2)
{
  c1 = na.omit(data1[,i])
  c2 = na.omit(data2[,i])
  rank_gen[i] = mean(c1)
  rank_epi[i] = mean(c2)
  col[i] = "darkgray"
  result = t.test(c1 ,c2,var.equal = FALSE)
  pval = result$p.value
  if(rank_gen[i] < rank_epi[i])
  {
    if(pval < 0.01)
    {
      col[i] = "red"
    }
  } else {
    if(pval < 0.005)
    {
      col[i] = "deepskyblue"
    }
  }
}


df <- data.frame(rank_gen, rank_epi, col)
df$col <- factor(df$col)

##I saved the df as a csv file, located in the Maintext_Data folder. Using that file - 
data = read.csv("Heterogeneity_outcomes_1B.csv") 

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
