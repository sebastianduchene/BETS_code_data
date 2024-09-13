setwd("~/Dropbox/projects_WORKING/BETS_prior_sensitivity/cholera_simulations_check/")

library(TreeSim)
library(NELSI)
library(phangorn)

for(i in 1:10){
  ntaxa <- 50
  # For the ultrametric tree
  tree <- sim.bd.taxa(n = ntaxa, lambda = 0.015, mu = 0.0025, numbsim = 1, complete = F)[[1]]
  tree
  # For a heterochronous tree use sim.bdsky.stt:
    #sim.bdsky.stt(n = 50, lambdasky = 0.06, deathsky = 0.02, sampprobsky = 0.2, timesky = 0)[[1]]
  tree <- read.tree(text = write.tree(tree))
  root_height <- max(allnode.times(tree))
  root_height
  
  sampling_times <- abs( round(sort( c(0,  rexp(ntaxa-1, 10 / root_height))) - 2009) )
  sampling_times
  
  #round(allnode.times(tree, tipsonly = T), 2)
  tree$tip.label <- paste(tree$tip.label, sampling_times, sep = "_")
  
  tree$root.edge 
  plot.tree.lines(tree, plot.new = T, line.type = "l", rotation.angle = 3*pi/2)
  
  # Simulate sampling times
  clock_rate <-  7e-7 * 2948589 / 5000
  clock_rate
  
  phylogram <- tree
  phylogram$edge.length <- phylogram$edge.length * clock_rate
  
  aln <- as.DNAbin(simSeq(phylogram, l = 5000))
  
  print(root_height)
  print(aln)
  print(length(seg.sites(aln)))
  print(max(sampling_times) - root_height)
  
  write.dna(aln, file = paste0("ultrametric_sample_cock_rate0.0004_sim", i, ".fasta"), format = "fasta", 
            nbcol = -1, colsep = "")
}

###################################
library(NELSI)
iqtree_tree <- read.tree("ultrametric_sample_cock_rate0.0004_sim10.fasta.treefile")
make.lsd.dates(iqtree_tree, random = F)

lsd_command <- "~/phyloApps/lsd-0.3beta/bin/lsd_unix -i ultrametric_sample_cock_rate0.0004_sim10.fasta.treefile -d out.date -r a -c"

system(lsd_command)
results_correct_rate <- gsub(".+ rate| , tMRCA.+", "",
                             grep("Tree.+ rate ", readLines("ultrametric_sample_cock_rate0.0004_sim10.fasta.treefile.result"), value = T))
results_correct_rate <- as.numeric(results_correct_rate)

randomised_rates <- vector()
for(i in 1:1000){
  make.lsd.dates(iqtree_tree, random = T)
  system(lsd_command)
  results_randomised_rate <- gsub(".+ rate| , tMRCA.+", "",
                               grep("Tree.+ rate ", readLines("ultrametric_sample_cock_rate0.0004_sim10.fasta.treefile.result"), value = T))
  randomised_rates[i] <- as.numeric(results_randomised_rate)
}

hist(randomised_rates, xlim = range(c(randomised_rates, results_correct_rate)), breaks = 50)
lines(rep(results_correct_date, 2), c(0, 500) , col = 'red', lwd = 5)
